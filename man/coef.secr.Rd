\name{coef.secr}
\alias{coef.secr}

\title{ Coefficients of secr Object }

\description{
  
  Extract coefficients (estimated beta parameters) from a spatially
  explicit capture--recapture model.
  
}

\usage{
\method{coef}{secr} (object, alpha = 0.05, ...)
}

\arguments{
  \item{object}{ \code{secr} object output from \code{secr.fit} }
  \item{alpha}{ alpha level }
  \item{...}{ other arguments (not used currently) }
}

\value{

A data frame with one row per beta parameter and columns for the
coefficient, SE(coefficient), asymptotic lower and upper 100(1--alpha)
confidence limits.

}

\seealso{ \code{\link{secr.fit}}, \code{\link{esaPlot}} }

\examples{

## load & extract coefficients of previously fitted null model
coef(secrdemo.0)

}

\keyword{ models }
