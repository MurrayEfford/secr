\name{capthist.parts}
\alias{animalID}
\alias{occasion}
\alias{trap}
\alias{alive}
\alias{xy}
\alias{alongtransect}
\alias{xy<-}
\alias{telemetryxy}
\alias{telemetryxy<-}
\alias{telemetered}

\title{ Dissect Spatial Capture History Object }

\description{
  Extract parts of an object of class `capthist'.
}

\usage{

animalID(object, names = TRUE, sortorder = c("snk", "ksn"))
occasion(object, sortorder = c("snk", "ksn"))
trap(object, names = TRUE, sortorder = c("snk", "ksn"))
alive(object, sortorder = c("snk", "ksn"))
alongtransect(object, tol = 0.01)
xy(object)
xy(object) <- value
telemetryxy(object, includeNULL = FALSE)
telemetryxy(object) <- value
telemetered(object)

}

\arguments{

  \item{object}{ a `capthist' object }
  \item{names}{ if FALSE the values returned are numeric indices rather
    than names }
    \item{sortorder}{character code for sort order (see Details)}
  \item{tol}{ tolerance for snapping to transect line (m) }
  \item{value}{ replacement value (see Details) }
  \item{includeNULL}{logical; if TRUE a NULL component is included for
  untelemetered animals } 

}
\details{

These functions extract data on detections, ignoring occasions when an
animal was not detected. By default, detections are ordered by occasion, animalID
and trap (\code{sortorder = "snk"}). The alternative is to order by 
trap, occasion and animalID (\code{sortorder = "ksn"}). (`n', `s' and `k' are the 
indices used internally for animals, occasions and traps respectively).

For historical reasons, "ksn" is used for locations within polygons and similar 
(\code{xy}). 

\code{trap} returns polygon or transect numbers if \code{traps(object)}
has detector type `polygon' or `transect'.

\code{alongtransect} returns the distance of each detection from the
start of the transect with which it is associated.

Replacement values must precisely match \code{object} in number of
detections and in their order. \code{xy<-} expects a dataframe of x and y
coordinates for points of detection within a `polygon' or `transect'
detector. \code{telemetryxy<-} expects a list of dataframes, one per telemetered animal.

}

\value{
For \code{animalID} and \code{trap} a vector of numeric or character values, one per detection.

For \code{alive} a vector of logical values, one per detection.

For \code{occasion}, a vector of numeric values, one per detection.

For \code{xy}, a dataframe with one row per detection and columns `x' and `y'.

If \code{object} has multiple sessions, the result is a list with one
component per session.

}

\seealso{ \code{\link{capthist}}, \code{\link{polyID}}, \code{\link{signalmatrix}} }

\examples{

## `captdata' is a demonstration dataset
animalID(captdata)

temp <- sim.capthist(popn = list(D = 1), make.grid(detector
    = "count"))
cbind(ID = as.numeric(animalID(temp)), occ = occasion(temp),
    trap = trap(temp))

}


