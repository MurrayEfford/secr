\name{sighting}
\alias{markocc}
\alias{markocc<-}
\alias{sighting}
\alias{sightings}
\alias{Tu}
\alias{Tu<-}
\alias{Tm}
\alias{Tm<-}
\alias{Tn}
\alias{Tn<-}

\title{ Sighting Attributes }

\description{
  Extract or replace the \code{markocc} attribute of a \code{traps}
  object that distinguishes marking occasions from sighting
  occasions. Also, extract or replace the attributes Tu, Tm and Tn of a capthist
  object, used for storing counts of sightings. All attributes are
  optional, but Tu, Tm and Tn require markocc to be specified.
}

\usage{
markocc(object, \dots)
markocc(object) <- value
sighting(object)
Tu(object, \dots)
Tu(object) <- value
Tm(object, \dots)
Tm(object) <- value
Tn(object, \dots)
Tn(object) <- value
}
\arguments{
  \item{object}{ \code{traps} object (\code{markocc}) or capthist object
  (\code{Tu}, \code{Tm}, \code{Tn})}
  \item{value}{ numeric matrix of detectors x occasions, or a vector
    (see Details) }
  \item{\dots}{ other arguments (not used) }
}
\details{
   
  For replacement of markocc, `value' should be a vector of integers
  indicating the occasions on which animals are sighted only (0) or
  marked or recaptured (1).

  For replacement of \code{Tu}, \code{Tm} or  \code{Tn}, `value' may be a scalar
  (total count) or a detectors x occasions matrix.  
  
}
\value{
  
  markocc(object) returns the markocc vector of the \code{traps}
  object. \code{markocc(object)} may be NULL.

  \code{Tu}, \code{Tm} and \code{Tn} return the respective attributes of a capthist object, or
  NULL if they are unspecified.

  sighting(object) returns TRUE if the markocc attribute indicates at least one 
  sighting-only occasion.
  
}

\seealso{
  
  \code{\link{traps}},
  \code{\link{addSightings}},
  \code{\link{sightingPlot}},
  \href{https://www.otago.ac.nz/density/pdfs/secr-markresight.pdf}{secr-markresight.pdf}
  
}


