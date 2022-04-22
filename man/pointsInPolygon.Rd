\name{pointsInPolygon}
\alias{pointsInPolygon}
\title{Points Inside Polygon}

\description{
Determines which of a set of  points lie inside a closed polygon or at
least one of a set of polygons
}

\usage{
pointsInPolygon(xy, poly, logical = TRUE)
}

\arguments{
  
  \item{xy}{2-column matrix or dataframe of x-y coordinates for points to assess}

  \item{poly}{2-column matrix or dataframe containing perimeter points of polygon, or
    a SpatialPolygonsDataFrame object from package \pkg{sp}, or a `mask'
    object (see Warning)}

  \item{logical}{logical to control the output when `poly' is a mask
  (see Details)}
  
}

\details{

  If \code{poly} is a SpatialPolygonsDataFrame object then the method
\code{over} is used from \pkg{sp}. This allows multiple polygons and
polygons with holes.

  If \code{poly} is an secr `mask' object then \code{xy} is discretized
  and matched to the cells in \code{poly}. If \code{logical = FALSE}
  then the returned value is a vector of integer indices to the row in
  `poly' corresponding to each row of `xy'; otherwise the result is a
  vector of logical values.
  
Otherwise, the algorithm is adapted from some code posted on the S-news list
by Peter Perkins (23/7/1996). The polygon should be closed (last point
same as first).

}

\value{
  Vector of logical or integer values, one for each row in \code{xy}
}

\seealso{
\code{\link[sp]{over}} 
}

\section{Warning}{If \code{poly} is a mask object then its cells must be
aligned to the x- and y- axes}

\examples{

## 100 random points in unit square
xy <- matrix(runif(200), ncol = 2)
## triangle centred on (0.5, 0.5)
poly <- data.frame(x = c(0.2,0.5,0.8,0.2), y = c(0.2,0.8,0.2,0.2))
plot(xy, pch = 1 + pointsInPolygon(xy, poly))
lines(poly)

}

\keyword{ manip }
