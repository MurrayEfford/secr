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
\alias{marked}
\alias{marked<-}

\title{ Sighting Attributes }

\description{
  Extract or replace the \code{markocc} attribute of a \code{traps}
  object that distinguishes marking occasions from sighting
  occasions. Also, extract or replace the attributes Tu, Tm and Tn of a capthist
  object, used for storing counts of sightings. All attributes are
  optional, but Tu, Tm and Tn require markocc to be specified. The attribute
  marked is the animal-level analogue of markocc: it records, for each animal and
  occasion, the fraction of the occasion for which the animal can be identified as
  marked, for example while a telemetry collar was working (telemetry type 'marking').
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
marked(object, \dots)
marked(object) <- value
}
\arguments{
  \item{object}{ \code{traps} object (\code{markocc}) or capthist object
  (\code{Tu}, \code{Tm}, \code{Tn}, \code{marked})}
  \item{value}{ numeric matrix of detectors x occasions, or a vector
    (see Details); for \code{marked} a matrix of animals x occasions }
  \item{\dots}{ other arguments (not used) }
}
\details{
   
  For replacement of markocc, `value' should be a vector of integers
  indicating the occasions on which animals are sighted only (0) or
  marked or recaptured (1).

  For replacement of \code{Tu}, \code{Tm} or  \code{Tn}, `value' may be a scalar
  (total count) or a detectors x occasions matrix.

  For replacement of \code{marked}, `value' is a numeric matrix with one row per
  animal and one column per occasion of the capthist object (telemetry occasions are
  ignored), or NULL. Each element is a known, non-negative multiplier of detector usage
  for that animal on that occasion, interpreted as the fraction of the occasion during
  which the animal is identifiable as marked: usually 0 (e.g. collar not working) or 1,
  or a fraction when the status changed during the occasion. The default (NULL) is 1 for
  all animals. In telemetry type 'marking' a collared animal contributes cues to the
  identified counts and Tm only while marked; at other times its cues are part of the
  unmarked sightings Tu. The values are known exposure, not probabilities: this is not
  availability in the sense of the \code{p.available} argument of \code{\link{sim.capthist}},
  which describes a random process unobserved in the data.

}
\value{
  
  markocc(object) returns the markocc vector of the \code{traps}
  object. \code{markocc(object)} may be NULL.

  \code{Tu}, \code{Tm}, \code{Tn} and \code{marked} return the respective attributes of a capthist object, or
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


