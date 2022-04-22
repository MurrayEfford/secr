\name{signal}
\alias{signalframe}
\alias{signal.capthist}
\alias{noise.capthist}
\alias{signalframe<-}
\alias{signal}
\alias{noise}
\alias{signal<-}
\alias{noise<-}

\title{ Signal Fields }

\description{
  
  Extract or replace signal attributes of a `capthist' object.
  
}

\usage{

signalframe(object)
signalframe(object) <- value

\method{signal}{capthist}(object, \dots)
\method{noise}{capthist}(object, \dots)

signal(object) <- value
noise(object) <- value
}

\arguments{

  \item{object}{ a `capthist' object }
  \item{value}{ replacement value (see Details) }
  \item{\dots}{ other arguments (not used)}
}
\details{

Signal attributes of a `capthist' object are stored in a dataframe
called the signalframe. This has one row per detection. The signalframe
includes the primary field `signal' and an unlimited number of other
fields. To extract the signal field alone use the signal method.

These functions extract data on detections, ignoring occasions when an
animal was not detected. Detections are ordered by occasion, animalID
and trap.

Replacement values must precisely match \code{object} in number of
detections and in their order.

}

\value{

For \code{signalframe} , a dataframe containing signal data and
covariates, one row per detection. The data frame has one row per
detection. See \code{\link{signalmatrix}} for a matrix with one row per
cue and columns for different microphones.

For \code{signal} and \code{noise}, a numeric vector with one element per detection.

If \code{object} has multiple sessions, the result is a list with one
component per session.

}

\seealso{ \code{\link{capthist}}, \code{\link{signalmatrix}} }

\examples{

## ovensong dataset has very simple signalframe
head(signalframe(signalCH))

}

\keyword{ manip }

