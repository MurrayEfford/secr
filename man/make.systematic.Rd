\name{make.systematic}
\alias{make.systematic}

\title{ Construct Systematic Detector Design }

\description{

A rectangular grid of clusters within a polygonal region.

}

\usage{

make.systematic(n, cluster, region, spacing = NULL, origin = NULL, 
    originoffset = c(0,0), chequerboard = c('all','black','white'), 
    order = c('x', 'y', 'xb', 'yb'), rotate = 0,  centrexy = NULL,
    keep.design = TRUE, ...)

}

\arguments{
  
  \item{n}{integer approximate number of clusters (see Details) }
  
  \item{cluster}{traps object defining a single cluster }
  
  \item{region}{dataframe or SpatialPolygonsDataFrame with coordinates
    of perimeter }

  \item{spacing}{ scalar distance between cluster centres }

  \item{origin}{vector giving x- and y-cooordinates of fixed grid origin
    (origin is otherwise random)}
    
  \item{originoffset}{ numeric; 2-vector (x,y offsets); see Details}
  
  \item{chequerboard}{ logical; if not `all' then alternate clusters are omitted }
  
  \item{order}{ character; sort order for clusters (see Details)}
  
  \item{rotate}{numeric; number of degrees by which to rotate entire design clockwise 
    about centroid of region bounding box}
    
  \item{centrexy}{numeric; 2-vector for centre of rotation, if any }   
  
  \item{keep.design}{logical; if TRUE then input argument values are retained}
   
  \item{\dots}{ arguments passed to \code{trap.builder} }
}

\details{

  \code{region} may be any shape. 
  
  \code{region} may be one of the spatial classes described 
  in \code{\link{boundarytoSF}}. Otherwise,
  \code{region} should be a dataframe with columns `x' and `y'.

  \code{spacing} may be a vector with separate values for spacing in x-
  and y- directions. If \code{spacing} is provided then \code{n} is ignored.
  
  If \code{n} is a scalar, the spacing of clusters is determined from
  the area of the bounding box of \code{region} divided by the requested
  number of clusters (this does not necessarily result in exactly n
  clusters). If \code{n} is a vector of two integers these are taken to be
  the number of columns and the number of rows.

  After preparing a frame of cluster centres, \code{make.systematic}
  calls \code{\link{trap.builder}} with method = `all'; \dots allows the
  arguments `rotation', `edgemethod', `plt', and `detector' to be
  passed. Setting the \code{trap.builder} arguments \code{frame},
  \code{method}, and \code{samplefactor} has no effect.
  
  Note the distinction between argument \code{rotate} and the \code{trap.builder} 
  argument \code{rotation} that is applied separately to each cluster.
  
  If \code{origin} is not specified then a random uniform origin is chosen within a box (width = spacing) placed with its bottom left corner at the bottom left of the bounding box of \code{region}, shifted by \code{originoffset}. Before version 3.1.8 the behaviour of \code{make.systematic} was equivalent to \code{originoffset = c(wx,wy)} where \code{wx,wy} are the cluster half widths.   
  
  \code{chequerboard = "black"} retains black `squares' and \code{chequerboard = "white"} retains white `squares', where the lower left cluster in the candidate rectangle of cluster origins is black, as on a chess board. The effect is the same as increasing spacing by sqrt(2) and rotating through 45 degrees.
  
  \code{order} determines the ordering of clusters in the resulting traps object. The options are a subset of those for \code{ID} argument of \code{\link{make.grid}}:
  
    \tabular{ll}{
     Option \tab Sort order \cr
     x \tab  column-dominant \cr
     y \tab  row-dominant \cr
     xb \tab  column-dominant boustrophedonical (alternate columns reversed) \cr
     yb \tab  row-dominant boustrophedonical (alternate rows reversed) \cr
     }

\code{rotate} rotates the array about the given centre (default is centroid of the bounding box of \code{region}). 

}

\value{

A single-session `traps' object.

From 3.2.0 these additional attributes are set --

\tabular{lll}{
origin \tab coordinates of grid origin \cr
centres \tab coordinates of true cluster centres (cf \code{\link{cluster.centres}}) \cr
originbox \tab vertices of rectangular spatial sampling frame for random origin \cr
}

From 4.2.0 if \code{keep.design} is TRUE then the input argument values are retained in attribute `design' (a list with first component \code{function = 'make.systematic'}).

}

\seealso{
  \code{\link{trap.builder}}, 
  \code{\link{make.lacework}},
  \code{\link{cluster.centres}}
  
}

\note{

    Do not confuse with the simpler function \code{\link{make.grid}},
    which places single detectors in a rectangular array.

}

\examples{

mini <- make.grid(nx = 2, ny = 2, spacing = 100)
region <- cbind(x=c(0,2000,2000,0), y=c(0,0,2000,2000))
temp <- make.systematic(25, mini, region, plt = TRUE)
temp <- make.systematic(c(6, 6), mini, region, plt = TRUE,
    rotation = -1)

## Example using shapefile "possumarea.shp" in
## "extdata" folder. By default, each cluster is 
## a single multi-catch detector

\dontrun{

library(sf)
shpfilename <- system.file("extdata/possumarea.shp", package = "secr")
possumarea <- st_read(shpfilename)

possumgrid <- make.systematic(spacing = 100, region =
    possumarea, plt = TRUE)

## or with 2 x 2 clusters
possumgrid2 <- make.systematic(spacing = 300,
    cluster = make.grid(nx = 2, ny = 2, spacing = 100),
    region = possumarea, plt = TRUE, edgemethod =
    "allinside")
## label clusters
text(cluster.centres(possumgrid2), levels(clusterID
    (possumgrid2)), cex=0.7)

## If you have GPSBabel installed and on the Path
## then coordinates can be projected and uploaded
## to a GPS with `writeGPS', which also requires the
## package `proj4'. Defaults are for a Garmin GPS
## connected by USB.

if (interactive()) {
    writeGPS(possumgrid, proj = "+proj=nzmg")
}

}

}


