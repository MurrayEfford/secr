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

The \dots argument is used to override previous settings (arguments of \code{\link{secr.fit}}); commonly method = "none" to recompute variances. A new 'details' argument may be incomplete: specified components replace those in the original (see example resetting 'saveprogress' below). 

The 'start' argument may not be overridden: it is set by \code{secr.refit} from the previous fit. The 'binomN' argument is elevated automatically from details$binomN if required.

A third use of \code{secr.refit} is to vary the model structure. This requires that 'object' is a fitted model. A change in model structure is implied when the \dots argument replaces any of 'model', 'CL', 'detectfn', 'link', 'fixed', 'hcov', or 'groups'. New starting values are then determined in \code{secr.fit} by a call to \code{\link{makeStart}}, and may not precisely follow the previous fit.

Also, \code{secr.refit} will compute a full model from a relative density model (i.e. estimate the intercept) - just set CL = FALSE. If method = "none" then new variances are computed using fdHess. Note this does \emph{not} work if the model includes individual covariates, which are prohibited in a full-likelihood model.

}

\value{

A fitted secr model.

}

\seealso{

\code{\link{secr.fit}},
\link{saveprogress}

}

\note{

In the case of a relative density model (i.e. a model was specified for D 
when CL = TRUE), a warning will be generated "overriding provided fixedbeta[1] for D". 
This can be ignored.

}

\examples{

\donttest{

# Repeat variance step using nlme::fdHess
# method = "none" skips likelihood maximization
secr.refit(secrdemo.0, method = "none")

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

# Example using secr.refit() to complete a relative density model
# fitCL <- secr.fit(OVpossumCH[[1]], buffer = 120, detectfn = 'HEX', 
#    CL = TRUE, model = D~x+y)
# fit <- secr.refit(fitCL, CL = FALSE, method = "none")    
# coef(fitCL)    
# coef(fit)
# predict(fit)

}

}
