\name{updateCH}
\alias{updateCH}

\title{Update Old capthist Format}

\description{

Before version 3.0, the internal data format for data from exclusive detectors (single, multi, proximityX, transectX) was a matrix with one row per detected animal and one column per sampling occasion; each cell was either zero or the number of the detector at which the animal was detected (with switched sign if the animal died). The format for data from proximity and other detectors was a 3-dimensional array (third dimension corresponding to detectors) that allowed more than one detection per animal per occasion.

From secr 3.0 all capthist data use the 3-D format internally. This simplifies a lot of the coding, and enables mixing of detector types within a session. The constraint that only one detection is allowed per animal per occasion at exclusive detectors is imposed by verify().

The data input functions (read.capthist etc.) automatically generate objects in the new format. Objects created and saved under earlier versions should be converted if they relate to the `exclusive' detector types listed above.

}

\usage{
updateCH(object)
}

\arguments{
  \item{object}{capthist object}
}

\details{
The function reduce.capthist is applied with the nominal detector type as the outputdetector.
}

\value{

Object with same class as the input. 

Updating has the side effect of discarding invalid supernumerary detections (e.g. if there were two detections of an animal on one occasion, only one will be included).

}

\examples{

# if we had the old ovenCH !
sapply(ovenCH, dim)
sapply(updateCH(ovenCH), dim)

}

