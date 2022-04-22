\name{rectangularMask}
\alias{rectangularMask}

\title{Rectangular Mask} 

\description{
Convert a mask or Dsurface with an irregular outline into a mask or
Dsurface with a rectangular outline and the same bounding box. This
enables contour plotting.
}

\usage{
rectangularMask(mask)
}

\arguments{
  \item{mask}{object of class mask or Dsurface}
}

\details{

The covariates of new points are set to missing. The operation may be
reversed (nearly) with \code{subset(rectmask, attr(rectmask, "OK"))}.

The results are unpredictable if the mask has been rotated.
}

\value{

  A rectangular mask or Dsurface with the same `area', `boundingbox',
  `meanSD', `polygon' and `polygon.habitat' attributes as \code{mask}. A
  logical vector attribute `OK' is added identifying the points
  inherited from \code{mask}.

}

\seealso{

  \code{\link{plot.Dsurface}}
  
}

\examples{

rMask <- rectangularMask(possummask)
plot(rMask)
plot(possummask, add = TRUE, col = "blue")

}
\keyword{ manip }

