\name{gridCells}
\alias{gridCells}

\title{Construct Grid Cells}

\description{

Forms grid cells centred on input points. 

}

\usage{

gridCells(x, cellsize = spacing(x), crs = NA)

}

\arguments{
  \item{x}{matrix or dataframe with x- and y-coordinates}
  \item{cellsize}{length of gridcell side}
  \item{crs}{crs description suitable for \code{\link{st_crs}} }
}

\details{

The argument x will often be a traps or mask object with spacing attribute. Otherwise \code{cellsize} must be provided.

}

\value{

A simple features (sf) object of class `sfc_MULTIPOLYGON'.

\code{crs} may be the integer  EPSG code (e.g. 3578 Yukon Albers).

}

\seealso{

\code{\link{plotMaskEdge}}, \code{\link{spacing}}

}

\examples{

plot(gridCells(traps(captdata)))
plot(traps(captdata), add = TRUE)
}

