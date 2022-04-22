\name{make.lacework}
\alias{make.lacework}

\title{Construct Lacework Detector Design}

\description{

A lacework design comprises a square grid with detectors placed at regular distances along the grid lines (Efford unpubl.). This requires fewer detectors than uniform coverage at close spacing and is simpler than clustered designs, while providing good spatial coverage and protection from alignment bias (Efford 2019).

}

\usage{

make.lacework(region, spacing = c(100, 20), times = NULL, origin = NULL, 
    rotate = 0, radius = NULL, detector = "multi", keep.design = TRUE)

}

\arguments{

\item{region}{dataframe or SpatialPolygonsDataFrame with coordinates of perimeter}

\item{spacing}{numeric 2-vector with major (grid) and minor spacings, or minor spacing only}

\item{times}{numeric ratio major:minor spacing if spacing length 1}
  
\item{origin}{numeric vector giving x- and y-cooordinates of fixed grid origin
    (origin is otherwise random)}

\item{rotate}{numeric; number of degrees by which to rotate design clockwise 
    about centroid of region bounding box}
    
\item{radius}{numeric; detectors are dropped if they are further than this 
    from a crossing}    

\item{detector}{character detector type -- see \code{\link{make.grid}}}   

\item{keep.design}{logical; if TRUE then input argument values are retained}

}

\details{

It is tidy for the major spacing (\code{spacing[1]}) to be a multiple of the minor spacing (\code{spacing[2]}); precisely one detector is then placed at each grid intersection. This outcome may also be achieved by providing only the minor spacing in the \code{spacing} argument and specifying an integer value for \code{times}.

In general it is better not to specify \code{origin}. Specifying both \code{origin} and \code{rotate} may result in incomplete coverage, as the desired grid is relative to the bounding box of the \emph{rotated} region.

Set \code{radius} < spacing[1]/2 to break lacework into multiple cross-shaped arrays centred on the intersections (crossing points) and truncated at \code{radius} metres (assuming you follow advice and express all linear measurements in metres).

}

\value{

An \pkg{secr} traps object. The attribute `crossings' is a 2-column matrix with the coordinates of the intersection points. If \code{keep.design} is TRUE then the input argument values are retained in attribute `design' (a list with first component \code{function = 'make.lacework'}).

}

\seealso{

\code{\link{make.systematic}}

}

\references{

Efford, M. G. (2019) Non-circular home ranges and the estimation of population density. \emph{Ecology} \bold{100}, e02580. https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.2580

}

\examples{

trps <- make.lacework(possumarea, c(1000,100), rotate = 45, detector = 'proximity')
plot(trps, gridspace = 1000)
lines(possumarea)
points(attr(trps, 'crossings'), pch = 16)

}

\keyword{ manip }