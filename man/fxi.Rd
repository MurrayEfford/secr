\name{fxi}
\alias{fxi.contour}
\alias{fxi.mode}
\alias{fxi.secr}
\alias{fxi}

\title{ Probability Density of Home Range Centre }

\description{
  Display contours of the probability density function for the estimated
  location of one or more range centres, compute values for
  particular points X, or compute mode of pdf. The pdf is given by
  \eqn{f(X_j|\omega_i) = \mbox{Pr}(\omega_i|X_j)\pi(X_j)}{f(X_j|\omega_i) =
  Pr(\omega_i|X)\pi(X_j)}, where \eqn{\pi(X)} is the probability density
  of range centres across the mask (Borchers and Efford 2008).
}

\usage{
fxi.contour (object, i = 1, sessnum = 1, border = 100, nx = 64,
    levels = NULL, p = seq(0.1,0.9,0.1), plt = TRUE, add = FALSE,
    fitmode = FALSE, plotmode = FALSE, fill = NULL,
    output = c('list','sf','SPDF'), ncores = NULL, ...)
fxi.secr(object, i = NULL, sessnum = 1, X = NULL, ncores = NULL)
fxi.mode(object, i = 1, sessnum = 1, start = NULL, ncores = NULL, ...)
}

\arguments{
  \item{object}{ a fitted secr model }
  \item{i}{ integer or character vector of individuals (defaults to all in fxi.secr),
  or a single individual as input to fxi.mode }
  \item{sessnum}{ session number if \code{object$capthist} spans
  multiple sessions}
  \item{border}{ width of blank margin around the outermost detectors }
  \item{nx}{ dimension of interpolation grid in x-direction }
  \item{levels}{ numeric vector of confidence levels for Pr(X|wi)}
  \item{p}{ numeric vector of contour levels as probabilities }
  \item{plt}{logical to plot contours}
  \item{add}{logical to add contour(s) to an existing plot }
  \item{fitmode}{logical to refine estimate of mode of each pdf}
  \item{plotmode}{logical to plot mode of each pdf}
  \item{X}{2-column matrix of x- and y- coordinates (defaults to mask)}
  \item{fill}{ vector of colours to fill contours (optional)}
  \item{output}{ character; format of output (list, sf or SpatialPolygonsDataFrame)}
  \item{ncores}{ integer number of threadss to be used for parallel processing}
  \item{start}{ vector of x-y coordinates for maximization }
  \item{\dots}{additional arguments passed to \code{contour} or \code{nlm}}
}

\details{

  \code{fxi.contour} computes contours of probability density for one
  or more detection histories. Increase \code{nx} for smoother
  contours. If \code{levels} is not set, contour levels are set
  to approximate the confidence levels in \code{p}.
  
  \code{fxi.secr} computes the probability density for one or more
  detection histories; \code{X} may contain coordinates for one or
  several points; a dataframe or vector (x then y) will be coerced to a
  matrix.

  \code{fxi.mode} attempts to find the x- and y-coordinates
  corresponding to the maximum of the pdf for a single detection history
  (i.e. \code{i} is of length 1). \code{fxi.mode} calls
  \code{\link[stats]{nlm}}.

  \code{fxi.contour} with \code{fitmode = TRUE} calls \code{fxi.mode}
  for each individual. Otherwise, the reported mode is an approximation
  (mean of coordinates of highest contour).

  If \code{i} is character it will be matched to row names of
  object$capthist (restricted to the relevant session in the case of a
  multi-session fit); otherwise it will be interpreted as a row number.
  
  Values of the pdf are normalised by dividing by the
  integral of \eqn{\mbox{Pr}(\omega_i|X)\pi(X)}{Pr(\omega_i|X)\pi(X)}
  over the habitat mask in \code{object}. (May differ in secr 4.0).
  
  Setting \code{ncores = NULL} uses the existing value from the environment variable 
  RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

  If \code{start} is not provided to \code{fit.mode} then (from 2.9.4) the weighted mean of 
  all detector sites is used (see Warning below).

  The \dots argument gives additional control over a contour plot; for
  example, set \code{drawlabels = FALSE} to suppress contour labels.
}

\value{
  
\code{fxi.contour} (output = 'list') --
  
Coordinates of the plotted contours are returned as a list with one
component per polygon. The list is returned invisibly if plt = TRUE.

An additional component `mode' reports the x-y coordinates of the
highest point of each pdf (see Details).

\code{fxi.contour} (output = 'SPDF') --
  
Contours are returned as a SpatialPolygonsDataFrame (see package
\pkg{sp}) with one component per animal. The attributes dataframe has two columns, the x-
and y-coordinates of the mode. The SpatialPolygonsDataFrame is returned
invisibly if plt = TRUE.

\code{fxi.contour} (output = 'sf') -- simple features 'sf' object, as for SPDF.

\code{fxi.secr} --

Vector of probability densities

\code{fxi.mode} --

List with components `x' and `y'

}
\references{

  Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

}

\seealso{
\code{\link{pdot.contour}}, \code{\link[graphics]{contour}}, \code{\link{fx.total}}
}

\note{

  From \pkg{secr} 2.8.3, these functions work with both homogeneous
  and inhomogeneous Poisson density models, and \code{fxi.secr} accepts
  vector-valued \code{i}.

  See \code{\link{fx.total}} for a surface summed across individuals.

}

\section{Warnings}{

  \code{fxi.mode} may fail to find the true mode unless a good starting
  point is provided. Note that the distribution may have multiple modes and 
  only one is reported. The default value of \code{start} before \pkg{secr} 2.9.4 
  was the first detected location of the animal. 

}

\examples{

\dontrun{

fxi.secr(secrdemo.0, i = 1, X = c(365,605))

## contour first 5 detection histories
plot(secrdemo.0$capthist)
fxi.contour (secrdemo.0, i = 1:5, add = TRUE,
    plotmode = TRUE, drawlabels = FALSE)

## extract modes only
## these are more reliable than those from fit.mode called directly as
## they use a contour-based approximation for the starting point
fxiout <- fxi.contour (secrdemo.0, i = 1:5, plt = FALSE, fitmode = TRUE)
t(sapply(fxiout, "[[", "mode"))

## using fill colours
## lty = 0 suppresses contour lines
## nx = 256 ensures smooth outline
plot(traps(captdata))
fxi.contour(secrdemo.0, i = 1:5, add = TRUE, p = c(0.5,0.95), drawlabels
    = FALSE, nx = 256, fill = topo.colors(4), lty = 0)

## output as SpatialPolygonsDataFrame
spdf <- fxi.contour(secrdemo.0, i = 1:3, plt = FALSE, p = c(0.5,0.95),
    nx = 256, output = 'SPDF', fitmode = TRUE)

## output as simple features
sf <- fxi.contour(secrdemo.0, i = 1:3, plt = FALSE, p = c(0.5,0.95),
    nx = 256, output = 'sf', fitmode = TRUE)

## save as ESRI shapefile test.shp etc.
library(sf)
st_write(sf, 'testsf.shp')

## plot contours and modes
plot(spdf)
points(data.frame(spdf))

}

}

\keyword{ hplot }