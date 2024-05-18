\name{traps.info}
\alias{searcharea}
\alias{polyID}
\alias{polyID<-}
\alias{transectID}
\alias{transectID<-}
\alias{transectlength}

\title{ Detector Attributes }

\description{
  Extract or replace attributes of an object of class `traps'.
}

\usage{
polyID(object)
polyID(object) <- value
transectID(object)
transectID(object) <- value
searcharea(object)
transectlength(object)
}

\arguments{
  \item{object}{ a `traps' object }
  \item{value}{ replacement value (see Details) }
}

\details{

The `polyID' and `transectID' functions assign and extract the attribute of a `traps'
object that relates vertices (rows) to particular polygons or transects. The
replacement value should be a factor of length equal to nrow(object).

The `searcharea' of a `polygon' traps object is a vector of the areas of
the component polygons in hectares. This value is read-only.

The `transectlength' of a `transect' traps object is a vector of the
lengths of the component transects in metres. This value is read-only.

}

\value{

\code{polyID} - a factor with one level per polygon. If the object does not have a polyID attribute then a factor with one level for each detector.
\code{searcharea} - numeric value of polygon areas, in hectares.
\code{transectlength} - numeric value of transect lengths, in metres.

}

\seealso{ 
  \code{\link{traps}}
}

\examples{

## default is a single polygon
temp <- make.grid(detector = "polygon", hollow = TRUE)
polyID(temp)
plot(temp)

## split in two
temp <- make.grid(detector = "polygon", hollow = TRUE)
polyID(temp) <- factor(rep(c(1,2),rep(10,2)))
plot(temp)
 
}

\keyword{ manip }

