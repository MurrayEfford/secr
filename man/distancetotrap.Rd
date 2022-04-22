\name{distancetotrap}
\alias{nearesttrap}
\alias{distancetotrap}

\title{ Distance To Nearest Detector }

\description{

Compute Euclidean distance from each of a set of points to the nearest
detector in an array, or return the sequence number of the detector
nearest each point.

}

\usage{

distancetotrap(X, traps)

nearesttrap(X, traps)

}

\arguments{
  \item{X}{ coordinates }
  \item{traps}{ traps object or 2-column matrix of coordinates }
}
\details{

\code{distancetotrap} returns the distance from each point in X to the
nearest detector in \code{traps}. It may be used to restrict the points
on a habitat mask.

For traps objects with polygon detector type (polygon, polygonX), and for SpatialPolygons, 
the function sf::st_distance is used internally(from \pkg{secr} 4.5.2).

}

\value{
\code{distancetotrap} returns a vector of distances (assumed to be in metres). 

\code{nearesttrap} returns the index of the nearest trap. 

}

\seealso{ \code{\link{make.mask}} }

\examples{

  ## restrict a habitat mask to points within 70 m of traps  
  ## this is nearly equivalent to using make.mask with the 
  ## `trapbuffer' option
  temptrap <- make.grid()
  tempmask <- make.mask(temptrap)
  d <- distancetotrap(tempmask, temptrap)
  tempmask <- subset(tempmask, d < 70)

}


\keyword{ manip }

