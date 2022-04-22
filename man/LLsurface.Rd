\name{LLsurface}
\alias{LLsurface}
\alias{LLsurface.secr}

\title{Plot Likelihood Surface}

\description{

LLsurface is a generic function to calculate log likelihood over a grid of values of 
two coefficients (beta parameters) from a fitted model and optionally make an 
approximate contour plot of the log likelihood surface.

A method is provided for secr objects.

}
\usage{

LLsurface(object, ...)

\method{LLsurface}{secr}(object, betapar = c("g0", "sigma"), xval = NULL,
    yval = NULL, centre = NULL, realscale = TRUE, plot = TRUE,
    plotfitted = TRUE, ncores = NULL, ...)
}

\arguments{
  \item{object}{fitted model, \code{secr} object output from \code{secr.fit}}
  \item{betapar}{character vector giving the names of two beta parameters}
  \item{xval}{vector of numeric values for x-dimension of grid}
  \item{yval}{vector of numeric values for y-dimension of grid}
  \item{centre}{  vector of central values for all beta parameters}
  \item{realscale}{logical. If TRUE input and output of x and y is on
        the untransformed (inverse-link) scale.}
  \item{plot}{logical. If TRUE a contour plot is produced}
  \item{plotfitted}{logical. If TRUE the MLE from \code{object} is shown on the plot (+)}
  \item{ncores}{integer number of threads for parallel processing}
  \item{\dots}{other arguments passed to \code{\link{contour}}}
}

\details{

\code{centre} is set by default to the fitted values of the beta
parameters in \code{object}. This has the effect of holding parameters
other than those in \code{betapar} at their fitted values.

If \code{xval} or \code{yval} is not provided then 11 values are set at
equal spacing between 0.8 and 1.2 times the values in \code{centre} (on
the `real' scale if \code{realscale} = TRUE and on the `beta' scale
otherwise).

Contour plots may be customized by passing graphical parameters through
the \dots argument.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

}

\value{

A matrix of the log likelihood evaluated at each
grid point (rows x, columns y), invisibly if \code{plot = TRUE}. 
Failed evaluations return NA.

}

\note{

  \code{LLsurface} works for named `beta' parameters rather than
  `real' parameters. The default \code{realscale = TRUE} only works for
  beta parameters that share the name of the real parameter to which
  they relate i.e. the beta parameter for the base level of the real
  parameter. This is because link functions are defined for real
  parameters not beta parameters.

  The contours are approximate because they rely on
  interpolation. See Examples for a more reliable way to compare the
  likelihood at the MLE with nearby points on the surface.
}

\examples{

\dontrun{

LLsurface(secrdemo.CL, xval = seq(0.16,0.40,0.02),
    yval = 25:35, nlevels = 20)

## now verify MLE
## click on MLE and apparent `peak'
if (interactive()) {
    xy <- locator(2)
    LLsurface(secrdemo.CL, xval = xy$x, yval = xy$y, plot = FALSE)
}

}
}

\keyword{ hplot }
