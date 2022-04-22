\name{contour}
\alias{pdot.contour}
\alias{buffer.contour}

\title{
Contour Detection Probability
}
\description{

  Display contours of the net probability of detection p.(X), or the
area within a specified distance of detectors.  \code{buffer.contour}
adds a conventional `boundary strip' to a detector (trap) array, where
\code{buffer} equals the strip width.

}
\usage{

pdot.contour(traps, border = NULL, nx = 64, detectfn = 0,
    detectpar = list(g0 = 0.2, sigma = 25, z = 1), noccasions = NULL,
    binomN = NULL, levels = seq(0.1, 0.9, 0.1), poly =
    NULL, poly.habitat = TRUE, plt = TRUE, add = FALSE, fill = NULL, ...)

buffer.contour(traps, buffer, nx = 64, convex = FALSE, ntheta = 100,
     plt = TRUE, add = FALSE, poly = NULL, poly.habitat = TRUE, 
     fill = NULL, ...)

}

\arguments{
  \item{traps}{ \code{traps} object (or mask for \code{buffer.contour}) }
  \item{border}{ width of blank margin around the outermost detectors }
  \item{nx}{ dimension of interpolation grid in x-direction }
  \item{detectfn}{ integer code or character string for shape of
    detection function 0 = halfnormal etc. -- see \link{detectfn}} 
  \item{detectpar}{ list of values for named parameters of detection
  function }
  \item{noccasions}{ number of sampling occasions}
  \item{binomN}{ integer code for discrete distribution (see
    \code{\link{secr.fit}}) }
  \item{levels}{ vector of levels for p.(X)}
  \item{poly}{ matrix of two columns, the x and y
    coordinates of a bounding polygon (optional) }
  \item{poly.habitat}{ logical as in \code{\link{make.mask}} }
  \item{plt}{ logical to plot contours }
  \item{add}{ logical to add contour(s) to an existing plot }
  \item{fill}{ vector of colours to fill contours (optional)}
  \item{\dots}{ other arguments to pass to \code{contour} }
  \item{buffer}{ vector of buffer widths }
  \item{convex}{ logical, if TRUE the plotted contour(s) will be convex }
  \item{ntheta}{ integer value for smoothness of convex contours }
}
\details{
\code{pdot.contour} constructs a rectangular mask and applies \code{\link{pdot}} to
compute the p.(X) at each mask point.

If \code{convex = FALSE}, \code{buffer.contour} constructs a mask and
contours the points on the basis of distance to the nearest detector at the
levels given in \code{buffer}.

If \code{convex = TRUE}, \code{buffer.contour} constructs a set of
potential vertices by adding points on a circle of radius =
\code{buffer} to each detector location; the desired contour is the
convex hull of these points (this algorithm derives from Efford, 2012).

If \code{traps} has a \link{usage} attribute then \code{noccasions} is
  set accordingly; otherwise it must be provided.

If \code{traps} is for multiple sessions then detectpar should be a list
of the same length, one component per session, and noccasions may be a
numeric vector of the same length.

Increase \code{nx} for smoother lines, at the expense of speed.

}
\value{
  
  Coordinates of the plotted contours are returned as a list with one
  component per polygon. The list is returned invisibly if \code{plt =
    TRUE}.

  For multi-session input (\code{traps}) the value is a list of such
  lists, one per session.
  
}

\references{

  Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially
explicit capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand
\url{https://www.otago.ac.nz/density/}.

}

\note{
 
  The precision (smoothness) of the fitted line in \code{buffer.contour}
  is controlled by \code{ntheta} rather than \code{nx} when \code{convex
    = TRUE}.

  To suppress contour labels, include the argument \code{drawlabels =
  FALSE} (this will be passed via \dots to \code{contour}). Other useful
  arguments of \code{contour} are \code{col} (colour of contour lines)
  and \code{lwd} (line width).

  You may wish to consider function st_buffer in package \pkg{sf} as an
  alternative to \code{buffer.contour}..

  \code{buffer.contour} failed with multi-session \code{traps} before
  \pkg{secr} 2.8.0.
  
}

\seealso{

  \code{\link{pdot}}, \code{\link{make.mask}}
  
}
\examples{

possumtraps <- traps(possumCH)

## convex and concave buffers
plot(possumtraps, border = 270)
buffer.contour(possumtraps, buffer = 100, add = TRUE, col = "blue")
buffer.contour(possumtraps, buffer = 100, convex = TRUE, add = TRUE)

## areas
buff.concave <- buffer.contour(possumtraps, buffer = 100,
    plt = FALSE)
buff.convex  <- buffer.contour(possumtraps, buffer = 100,
    plt = FALSE, convex = TRUE)
sum (sapply(buff.concave, polyarea)) ## sum over parts
sapply(buff.convex, polyarea)

## effect of nx on area
buff.concave2 <- buffer.contour(possumtraps, buffer = 100,
    nx = 128, plt = FALSE)
sum (sapply(buff.concave2, polyarea))

\dontrun{

plot(possumtraps, border = 270)
pdot.contour(possumtraps, detectfn = 0, nx = 128, detectpar =
    detectpar(possum.model.0), levels = c(0.1, 0.01, 0.001),
    noccasions = 5, add = TRUE)

## clipping to polygon
olddir <- setwd(system.file("extdata", package = "secr"))
possumtraps <- traps(possumCH)
possumarea <- read.table("possumarea.txt", header = TRUE)
par(xpd = TRUE, mar = c(1,6,6,6))
plot(possumtraps, border = 400, gridlines = FALSE)
pdot.contour(possumtraps, detectfn = 0, nx = 256, detectpar =
    detectpar(possum.model.0), levels = c(0.1, 0.01, 0.001),
    noccasions = 5, add = TRUE, poly = possumarea, col = "blue")
lines(possumarea)
setwd(olddir)
par(xpd = FALSE, mar = c(5,4,4,2) + 0.1)    ## reset to default

}

}

\keyword{ hplot }
