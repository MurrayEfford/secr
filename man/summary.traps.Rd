\name{summary.traps}
\alias{summary.traps}
\alias{print.summary.traps}

\title{ Summarise Detector Array }

\description{
  Concise description of \code{traps} object.
}

\usage{
\method{summary}{traps}(object, getspacing = TRUE, covariates = FALSE, \dots)
\method{print}{summary.traps}(x, terse = FALSE, \dots)
}

\arguments{
  \item{object}{    \code{\link{traps}} object }
  \item{getspacing}{ logical to calculate spacing of detectors from scratch }
  \item{covariates}{ logical; if true each covariate is summarised}
  \item{x}{ \code{summary.traps} object }
  \item{terse}{ if TRUE suppress printing of usage and covariate summary }
  \item{\dots}{ arguments passed to other functions }
}

\value{
An object of class \code{summary.traps}, a list with elements
  \item{detector }{\code{\link{detector}} type ("single", "multi", "proximity" etc.)}
  \item{ndetector }{number of detectors}
  \item{xrange }{range of x coordinates }
  \item{yrange }{range of y coordinates }
  \item{spacing }{mean distance from each trap to nearest other trap }
  \item{usage }{table of usage by occasion}
  \item{covar }{summary of covariates}
}

\details{

When \code{object} includes both categorical (factor) covariates and
\code{usage}, usage is tabulated for each level of the covariates. 

Computation of \code{spacing} (mean distance to nearest trap) is slow
and may hit a memory limit when there are many traps. In this case, turn
off the computation with \code{getspacing} = FALSE. 

}

\seealso{ \code{\link{print}}, \code{\link{traps}} }

\examples{
demo.traps <- make.grid()
summary(demo.traps)   ## uses print method for summary.traps object
}
\keyword{ models }

