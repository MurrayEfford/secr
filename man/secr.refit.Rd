\name{secr.refit}
\alias{secr.refit}

\title{Repeat or Resume Model Fit}

\description{

From 5.2.1 it is possible to save a part-fitted model to a file and resume 
fitting with \code{\link{secr.refit}}. The function may also be used to re-fit 
a completed model, or to repeat just the variance estimation (Hessian) step.

}

\usage{

secr.refit(object, ...)

}

\arguments{
  \item{object}{secr object (a fitted model) or character value naming a saved 
  progress RDS file (e.g. progress.RDS)}
  \item{\dots}{arguments of \code{\link{secr.fit}} to override those in object}
}

\details{

There are two common uses for the function.

1. To resume fitting of a partly fitted model from an RDS file saved with \code{\link{saveprogress}}.

2. To recompute variances of a fitted model. The numerical method is fdHess from package \pkg{nlme}. Often this delivers variance estimates when the default algorithm failed.

}

\value{

A fitted secr model.

}

\seealso{

\code{\link{secr.fit}},
\link{saveprogress}

}

\examples{

# Repeat variance step using nlme::fdHess
# method = "none" skips likelihood maximization
secr.refit(secrdemo.0, method = "none")

\donttest{

# Save data and coefficients after every 10 evaluations to a file
# with the default name 'progress.RDS'.
fit1 <- secr.fit(captdata, CL = TRUE, details = list(saveprogress = 10))

# Review progress. In this example, fit1 ran to completion, 
# but the log stops at a multiple of 10 evaluations.
attr(readRDS('progress.RDS'), 'log')

# Restart at last line in log, suppressing overwrite of progress log
secr.refit('progress.RDS', details = list(saveprogress = FALSE))

# cleanup
unlink("progress.RDS")

}

}
