\name{polyarea}
\alias{polyarea}

\title{
Area of Polygon(s)
}
\description{
Area of a single closed polygon (simple x-y coordinate input) or of
multiple polygons, possibly with holes.
}
\usage{
polyarea(xy, ha = TRUE)
}

\arguments{
  \item{xy}{dataframe or list with components `x' and `y', or a
    SpatialPolygons or SpatialPolygonsDataFrame object from package \pkg{sp}, or an \pkg{sf} object with polygon data}
  \item{ha}{logical if TRUE output is converted from square metres to hectares}
}
\details{
  For sf, sfc, SpatialPolygons or SpatialPolygonsDataFrame objects, the package \pkg{sf} is used.
}

\value{
  A scalar.
}

\seealso{
  \code{\link{bufferContour}}
}

\examples{
polyarea(make.grid(hollow = TRUE))
}

