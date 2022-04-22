\name{PG}
\alias{PG}

\title{Telemetry Fixes in Polygons}

\description{
  
For a telemetry dataset, either as a standalone capthist object with
detector type `telemetryonly' or the xylist attribute of a combined capthist
object resulting from \code{addTelemetry}, determine the proportion of
fixes of each individual that lie within a set of polygons. Typically
used to obtain the proportion of fixes on a trapping grid, hence
`proportion on grid'.

}

\usage{

PG(CH, poly = NULL, includeNULL = FALSE, plt = FALSE, ...)

}

\arguments{
  \item{CH}{ capthist object including telemetry locations}
  \item{poly}{ polygon object (see \code{\link{boundarytoSF}}) }
  \item{includeNULL}{ logical; if TRUE then missing values are returned
    for animals without telemetry data}
  \item{plt}{ logical; if TRUE then poly and telemetry locations are plotted}
  \item{\dots}{ other arguments passed to \code{\link{buffer.contour}}}
}

\details{
  
  By default \code{poly} is obtained by applying
  \code{\link{buffer.contour}} with arguments \dots to the traps
  attribute of \code{CH}. Note that either a positive \code{buffer}
  argument or \code{convex = TRUE} is needed for the polygon to have
  area > 0.

  If \code{plt = TRUE}, \code{buffer.contour} is used to plot
  \code{poly} and the points are overplotted (open circles outside,
  filled circles inside). To control the framing of the plot, create an
  initial plot (e.g., with plot.traps, setting the \code{border}
  argument) and use \code{add = TRUE} (see Examples).
  
}
\value{
  
  Numeric vector of proportions. If \code{includeNULL = TRUE} length
  equal to number of animals (rows) in \code{CH}; otherwise length is
  the number of animals for which there is telemetry data (because
  xylist may cover only a subset of animals in \code{CH}).

}

\seealso{
  
 \code{\link{addTelemetry}},  \code{\link{buffer.contour}}, 
\code{\link{pointsInPolygon}}

}

\references{
  
Grant, T. J. and Doherty, P. F. (2007) Monitoring of the flat-tailed
horned lizard with methods incorporating detection
probability. \emph{Journal of Wildlife Management} \bold{71}, 1050--1056

}


\examples{

\dontrun{
olddir <- setwd('d:/density communication/combining telemetry and secr/possums')
CvilleCH <- read.capthist('CVILLE summer captures 4occ.txt',
                          'CVILLE detectors summer 4occ.txt',
                          detector = 'single')
CvilleGPS <- read.telemetry('CVILLE GPS Combined 4occ.txt')
CvilleGPSnew <- read.telemetry('CVILLE summer GPS New occasions.txt')
setwd(olddir)

CvilleBoth <- addTelemetry(CvilleCH, CvilleGPSnew)
plot(CvilleBoth, border = 400)
PG(CvilleBoth, buffer = 100, convex = TRUE, plt = TRUE, add = TRUE, 
    col = 'red')

###################################################################
## this code computes an area-adjusted density estimate
## cf Grant and Doherty 2007
PGD <- function (CH, estimator = 'h2', ...) {
    pg <- PG(CH, ...)
    PGbar <- mean(pg)
    N <- closedN(CH, estimator)
    A <- polyarea(buffer.contour(traps(CH), ...)[[1]])
    Dhat <- N$Nhat / A * PGbar
    varDhat <- (N$Nhat^2 * var(pg) + PGbar^2 * N$seNhat^2) / A^2 
    c(Dhat = Dhat, seDhat = sqrt(varDhat))
}
plot(traps(CvilleBoth), border = 400)
PGD(CvilleBoth, buffer = 0, convex = TRUE, plt = TRUE, add = TRUE)
PGD(CvilleBoth, est='null', buffer = 0, convex = TRUE, plt = FALSE)

###################################################################
## this code generates a PG summary for telemetry records randomly
## translated and rotated, keeping the centres within a habitat mask

randomPG <- function(CH, poly = NULL, mask, reorient = TRUE, nrepl = 1,
                     seed = 12345, ...) {
    moveone <- function(xy, newcentre) {
        xy <- sweep(xy,2,apply(xy,2,mean))
        if (reorient)  ## random rotation about centre
            xy <- rotate(xy, runif(1)*360)
        sweep(xy,2,unlist(newcentre), "+")
    }
    onerepl <- function(r) {   ## r is dummy for replicate
        centres <- sim.popn(D = D, core = mask, model2D = "IHP",
                            Ndist = "fixed")
        xyl <- mapply(moveone, xyl, split(centres, rownames(centres)))
        attr(CH, 'xylist') <- xyl  ## substitute random placement
        PG(CH = CH , poly = poly, plt = FALSE, ...)
    }
    set.seed(seed)
    if (!requireNamespace('sf')) stop ("requires package sf")
    if (is.null(poly)) {
        poly <- buffer.contour (traps(CH), ...)
        poly <- lapply(poly, as.matrix)
        poly <- sf::st_sfc(sf::st_polygon(poly))
    }
    xyl <- telemetryxy(CH)
    D <- length(xyl) / maskarea(mask)
    sapply(1:nrepl, onerepl)
}

mask <- make.mask (traps(CvilleBoth), buffer = 400, type = "trapbuffer")
pg <- randomPG (CvilleBoth, mask = mask, buffer = 100, convex = TRUE,
    nrepl = 20)
apply(pg, 1, mean)
###################################################################

}
}

\keyword{ manip }

