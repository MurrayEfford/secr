\name{addCovariates}
\alias{addCovariates}

\title{
  Add Covariates to Mask or Traps
}

\description{

Tools to construct spatial covariates for existing mask or traps objects
from a spatial data source.  

}

\usage{

addCovariates(object, spatialdata, columns = NULL, strict = FALSE, replace = FALSE)

}

\arguments{
  \item{object}{mask, traps or popn object}
  \item{spatialdata}{spatial data source (see Details)}
  \item{columns}{character vector naming columns to include (all by default)}
  \item{strict}{logical; if TRUE a check is performed for points in \code{object} that
lie outside \code{spatialdata} (mask data sources only)}
  \item{replace}{logical; if TRUE then covariates with duplicate names are replaced; otherwise a new column is added}
}

\details{

  The goal is to obtain the value(s) of one or more spatial covariates
  for each point (i.e. row) in \code{object}. The procedure depends on
  the data source \code{spatialdata}, which may be either a spatial
  coverage (raster or polygon) or an object with covariate values at
  points (another mask or traps object). In the first case, an
  overlay operation is performed to find the pixel or
  polygon matching each point. In the second case, a search is conducted
  for the closest point in \code{spatialdata}.

  If \code{spatialdata} is a character value then it is interpreted as
  the name of a polygon shape file (excluding `.shp').

  If \code{spatialdata} is a SpatialPolygonsDataFrame, 
  SpatialGridDataFrame or 'sf' object from \pkg{sf} then it will be used 
  in an overlay operation as described. 

   If package \pkg{terra} has been installed then \code{spatialdata} may also be a RasterLayer from package \pkg{raster} or SpatRaster from \pkg{terra}. If provided \code{counts} should be a single name that will be used for the values (otherwise 'raster' will be used).
 
  If \code{spatialdata} is a \code{mask} or \code{traps} object then it
  is searched for the closest point to each point in \code{object}, and
  covariates are drawn from the corresponding rows in
  \code{covariates(spatialdata)}. By default (\code{strict = FALSE}),
  values are returned even when the points lie outside any cell of the mask.
  
}

\value{

  An object of the same class as \code{object} with new or augmented
  \code{covariates} attribute. Column names and types are derived from the input.

}

\seealso{

  \code{\link{make.mask}}, \code{\link{read.mask}}, \code{\link{read.traps}}

}

\section{Warning}{

  Use of a SpatialGridDataFrame for \code{spatialdata} is untested.
  
}

\examples{

## In the Lake Station skink study (see ?skink), habitat covariates were
## measured only at trap sites. Here we extrapolate to a mask, taking
## values for each mask point from the nearest trap.

LSmask <- make.mask(LStraps, buffer = 30, type = "trapbuffer")
tempmask <- addCovariates(LSmask, LStraps)
## show first few lines
head(covariates(tempmask))

}

