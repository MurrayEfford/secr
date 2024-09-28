\name{ellipse.secr}
\alias{ellipse.secr}
\alias{ellipse.bvn}

\title{ Confidence Ellipses }

\description{
  
Plot joint confidence ellipse for two parameters of secr model, or for a
distribution of points.

}

\usage{

ellipse.secr(object, par = c("g0", "sigma"), alpha = 0.05,
    npts = 100, plot = TRUE, linkscale = TRUE, add = FALSE,
    col = palette(), ...)

ellipse.bvn(xy, alpha = 0.05, npts = 100, centroid = TRUE, add = FALSE, ...)

}

\arguments{
  \item{object}{\code{secr} object output from \code{secr.fit}}
  \item{par}{character vector of length two, the names of two `beta' parameters}
  \item{alpha}{alpha level for confidence intervals }
  \item{npts}{number of points on perimeter of ellipse}
  \item{plot}{logical for whether ellipse should be plotted}
  \item{linkscale}{logical; if FALSE then coordinates will be backtransformed from the link scale}
  \item{add}{logical to add ellipse to an existing plot }
  \item{col}{vector of one or more plotting colours}
  \item{\dots}{arguments to pass to plot functions (or polygon() in the case
  of ellipse.bvn)}
  \item{xy}{2-column matrix of coordinates}
  \item{centroid}{logical; if TRUE the plotted ellipse is a confidence
    region for the centroid of points in \code{xy}}
}
\details{

\code{ellipse.secr} calculates coordinates of a confidence ellipse from
the asymptotic variance-covariance matrix of the beta parameters
(coefficients), and optionally plots it.

If \code{linkscale} == FALSE, the inverse of the appropriate link
transformation is applied to the coordinates of the ellipse, causing it
to deform.

If \code{object} is a list of secr models then one ellipse is
constructed for each model. Colours are recycled as needed.

\code{ellipse.bvn} plots a bivariate normal confidence ellipse for the
centroid of a 2-dimensional distribution of points (default centroid =
TRUE), or a Jennrich and Turner (1969) elliptical home-range model.

}
\value{
A list containing the x and y coordinates is returned invisibly from
either function.
}

\references{

  Jennrich, R. I. and Turner, F. B. (1969) Measurement of
   non-circular home range. \emph{Journal of Theoretical Biology},
  \bold{22}, 227--237.

}

\examples{

ellipse.secr(secrdemo.0)

}

