\name{raster}

\docType{methods}

\alias{raster}
\alias{raster,mask-method}
\alias{raster,Dsurface-method}

\alias{rast}
\alias{rast,mask-method}
\alias{rast,Dsurface-method}

\title{Create a RasterLayer Object from Mask or Dsurface}

\description{

  Methods to convert \pkg{secr} object to a RasterLayer object.

}


\usage{

\S4method{raster}{mask}(x, covariate, values = 1, crs = NA)

\S4method{raster}{Dsurface}(x, covariate, values = 1, crs = NA)
\S4method{rast}{mask}(x, covariate, values = 1, crs = "")
\S4method{rast}{Dsurface}(x, covariate, values = 1, crs = "")

}

\arguments{

  \item{x}{mask or Dsurface object}

  \item{covariate}{character name of covariate to provide values for
  RasterLayer}
  
  \item{values}{numeric values for RasterLayer}

  \item{crs}{character or object of class CRS. Optional PROJ.4 type description of
  a Coordinate Reference System (map projection). 
  
}

}

\details{

There are two ways to specify the values to be used. If \code{covariate} is provided then 
the values of the corresponding covariate of the mask or Dsurface are used. Otherwise, 
\code{values} is duplicated to the required number of rows.

The resulting RasterLayer may optionally include a PROJ.4 map projection
defined via \code{crs}. The specification may be very simple (as in the
example below) or complex, including an explicit datum and other
arguments. Projections are used by \pkg{sf}, \pkg{terra}, \pkg{raster}, 
\pkg{sp} and other packages. See \code{\link[raster]{raster}} for 
further explanation and links.

The S3 classes `mask' and `Dsurface' are defined in \pkg{secr} as
virtual S4 classes. This enables these extensions to the list of S4
methods defined in \pkg{raster} and \pkg{terra}.
  
Although these methods work `standalone', it is currently necessary to
load the \pkg{terra} or \pkg{raster} package to do much with the result (e.g., plot it).

}

\value{
  
RasterLayer (raster) \cr
SpatRaster (rast)

}

\seealso{
  
\code{\link[raster]{raster}}, \code{\link[terra]{rast}}

}

\note{

  Prior to \pkg{secr} 2.9.5 these methods could fail unpredictably
  because an intermediate array was badly dimensioned due to truncation of
  a floating point value.
  
}

\examples{

\dontrun{

shorePossums <- predictDsurface(possum.model.Ds)
tmp <- raster(shorePossums, covariate = "D.0")
library(raster)
plot(tmp, useRaster = FALSE)

## alternative with same result
tmp <- raster(shorePossums, values = covariates(shorePossums)$D.0)

## set the projection
## here the crs PROJ.4 spec refers simply to the old NZ metric grid
tmp <- raster(shorePossums, "D.0", crs = "+proj=nzmg")
## check the projection
proj4string(tmp)
     
# mask to rast
dtsrast <- rast(possummask, covariate = 'd.to.shore')
plot(dtsrast)

# Dsurface to rast
class(shorePossums) <- c('mask', 'data.frame') # or change "mask" to "Dsurface" above
Drast <- rast(shorePossums, covariate = 'D.0')
plot(Drast)

}

}

