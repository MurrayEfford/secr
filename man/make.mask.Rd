\name{make.mask}
\alias{make.mask}

\title{ Build Habitat Mask }

\description{
Construct a habitat mask object for spatially explicit capture-recapture. 
A mask object is a set of points with optional attributes. 
}

\usage{

make.mask(traps, buffer = 100, spacing = NULL, nx = 64, ny = 64, type =
    c("traprect", "trapbuffer", "pdot", "polygon", "clusterrect",
    "clusterbuffer", "rectangular", "polybuffer"), poly = NULL, poly.habitat = TRUE,
    cell.overlap = c("centre","any","all"), keep.poly = TRUE, check.poly = TRUE, 
    pdotmin = 0.001, random.origin = FALSE,
    ...)

}

\arguments{
  \item{traps}{ object of class \code{traps} }
  \item{buffer}{ width of buffer in metres }
  \item{spacing}{ spacing between grid points (metres) }
  \item{nx}{ number of grid points in `x' direction }
  \item{ny}{ number of grid points in `y' direction (type = `rectangular') }
  \item{type}{ character string for method }
  \item{poly}{ bounding polygon to which mask should be clipped (see Details) }
  \item{poly.habitat}{ logical for whether poly represents habitat or
  its inverse (non-habitat)}
  \item{cell.overlap}{character string for cell vertices used to determine overlap with `poly'}
  \item{keep.poly}{ logical; if TRUE any bounding polygon is saved as the
  attribute `polygon' }
  \item{check.poly}{ logical; if TRUE a warning is given for traps that
  lie outside a bounding polygon }
  \item{pdotmin}{ minimum detection probability for inclusion in mask when \code{type = "pdot"} (optional) }
  \item{random.origin}{ logical; if TRUE the mask coordinates are jittered}
  \item{...}{ additional arguments passed to \code{pdot} when \code{type = "pdot"} }
}

\details{

The `traprect' method constructs a grid of points in the rectangle
formed by adding a buffer strip to the minimum and maximum x-y
coordinates of the detectors in \code{traps}. Both `trapbuffer' and
`pdot' start with a `traprect' mask and drop some points.

The `trapbuffer' method restricts the grid to points within distance
\code{buffer} of any detector.

The `pdot' method restricts the grid to points for which the net
detection probability \eqn{p.(\mathbf{X})}{p.(X)} (see
\code{\link{pdot}}) is at least \code{pdotmin}. Additional parameters
are used by \code{pdot} (detectpar, noccasions). Set these with the
\dots argument; otherwise \code{make.mask} will silently use the
arbitrary defaults. \code{pdot} is currently limited to a halfnormal
detection function.

The `clusterrect' method constructs a grid of rectangular submasks
centred on `clusters' of detectors generated with
\code{\link{trap.builder}} (possibly indirectly by
\code{\link{make.systematic}}). The `clusterbuffer' method resembles
`trapbuffer', but is usually faster when traps are arranged in clusters
because it starts with a `clusterrect' mask.

The `rectangular' method constructs a simple rectangular mask with the
given nx, ny and spacing. 

The `polybuffer' method constructs a mask by buffering around the polygon specified in the `poly' argument. If that inherits from `SpatialPolygons' or `sfc' then the buffering is performed with sf::st_buffer. Otherwise, buffering is approximate, based on the distance to points on an initial discretized mask enclosed by `poly' (points at half the current `spacing').

If \code{poly} is specified, points outside \code{poly} are dropped (unless type = "polybuffer"). The default is to require only the centre to lie within \code{poly}; use \code{cell.overlap = "all"} to require all cell corners to lie within \code{poly}, or \code{cell.overlap = "any"} to accept cells with any corner in \code{poly}. The `polygon' method places points on a rectangular grid clipped to the
polygon (\code{buffer} is not used). Thus `traprect' is equivalent to
`polygon' when \code{poly} is supplied. \code{poly} may be either
\itemize{
\item a matrix or dataframe of two columns interpreted as x and y
coordinates, or
\item an object from package `sf' with polygon geometries, or
\item a SpatialPolygons or SpatialPolygonsDataFrame object as defined in the package `sp',
possibly imported by reading a shapefile.
}

If \code{spacing} is not specified then it is determined by dividing the
range of the x coordinates (including any buffer) by \code{nx}.

\code{random.origin} shifts the origin of the mask by a uniform random displacement within a spacing x spacing grid cell, while ensuring that the mask also satisfies the buffer requirement. \code{random.origin} is available only for `traprect', `trapbuffer', `polygon', and `rectangular' types, and \code{spacing} 
must be specified.

}

\value{
  
    An object of class \code{mask}.  When \code{keep.poly = TRUE},
    \code{poly} and \code{poly.habitat} are saved as attributes of the
    mask.

  }

\note{

A warning is displayed if \code{type = "pdot"} and the buffer is too small to
include all points with \eqn{p.} > pdotmin.

A habitat mask is needed to fit an SECR model and for some related
computations. The default mask settings in \code{secr.fit} may be good
enough, but it is preferable to use \code{make.mask} to construct a mask
in advance and to pass that mask as an argument to \code{secr.fit}.

The function \code{bufferContour} displays the extent of one or more
`trapbuffer' zones - i.e. the effect of buffering the detector array
with varying strip widths. 

}

\seealso{

  \code{\link{mask}}, \code{\link{read.mask}}, \code{\link{subset.mask}},
  \code{\link{pdot}}, \code{\link{bufferContour}}, \code{\link{deleteMaskPoints}},
  \code{\link{as.mask}}

}

\examples{

temptrap <- make.grid(nx = 10, ny = 10, spacing = 30)

## default method: traprect
tempmask <- make.mask(temptrap, spacing = 5)
plot(tempmask)
summary (tempmask)

## make irregular detector array by subsampling 
## form mask by `trapbuffer' method
temptrap <- subset (temptrap, sample(nrow(temptrap), size = 30))
tempmask <- make.mask (temptrap, spacing = 5, type = "trapbuffer")
plot (tempmask)
plot (temptrap, add = TRUE)

\dontrun{

## form mask by "pdot" method
temptrap <- make.grid(nx = 6, ny = 6)
tempmask <- make.mask (temptrap, buffer = 150, type = "pdot", 
    pdotmin = 0.0001, detectpar = list(g0 = 0.1, sigma = 30),
    noccasions = 4)
plot (tempmask)
plot (temptrap, add = TRUE)

## Using an ESRI polygon shapefile for clipping (shapefile
## polygons may include multiple islands and holes).

library(sf)
shpfilename <- system.file("extdata/possumarea.shp", package = "secr")
possumarea <- st_read(shpfilename)

possummask2 <- make.mask(traps(possumCH), spacing = 20,
    buffer = 250, type = "trapbuffer", poly = possumarea)
par(mar = c(1,6,6,6), xpd = TRUE)
plot (possummask2, ppoly = TRUE)
plot(traps(possumCH), add = TRUE)
par(mar = c(5,4,4,2) + 0.1, xpd = FALSE)

## if the polygon delineates non-habitat ...
seaPossumMask <- make.mask(traps(possumCH), buffer = 1000, 
    type = "traprect", poly = possumarea, poly.habitat = FALSE)
plot(seaPossumMask)
plot(traps(possumCH), add = TRUE)
## this mask is not useful!

}

%% ## OR, as this is a simple polygon, just...
%% possumarea <- possumarea@polygons[[1]]@Polygons[[1]]@coords
%% possummask2 <- make.mask(traps(possumCH), spacing = 20, poly = possumarea,
%%     buffer = 250, type = "trapbuffer")

}

