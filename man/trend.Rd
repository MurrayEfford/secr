\name{Trend}
\alias{predictDlambda}

\title{ Density Trend }

\description{

  Functions for multi-session density trend analysis.

}

\usage{

predictDlambda(object, alpha = 0.05)

}

\arguments{

  \item{object}{multi-session secr object output from secr.fit}
  \item{alpha}{alpha level for confidence intervals}

}

\details{
  
Usage is described in \href{https://www.otago.ac.nz/density/pdfs/secr-trend.pdf}{secr-trend.pdf}. Briefly, setting details argument 'Dlambda' causes the density model (D~xxx) to be interpreted as a session-specific trend model with parameters for the initial density (D1) and each subsequent session-on-session change in density \eqn{\lambda[t] = D[t+1]/D[t]}. 
  
}

\value{

A table of session-specific estimates (initial D, subsequent \eqn{\lambda[t]}) with SE and confidence intervals.

}

\seealso{
  
\code{\link{predictDsurface}}, 
\code{\link{secr.fit}}

}

\examples{

\donttest{
# a model with constant lambda
msk <- make.mask(traps(ovenCH[[1]]), buffer = 300, nx = 25)
fit <- secr.fit(ovenCH, model = D~1, mask = msk, trace = FALSE, 
                 details = list(Dlambda = TRUE), ncores = 2)
predictDlambda(fit)

}

}
\keyword{ classes }
