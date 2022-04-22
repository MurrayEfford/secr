\name{addTelemetry}
\alias{addTelemetry}
\alias{xy2CH}
\alias{telemetrytype}
\alias{telemetrytype<-}

\title{Combine Telemetry and Detection Data}

\description{
  
Animal locations determined by radiotelemetry can be used to augment
capture--recapture data.  The procedure in \pkg{secr} is first to form a
capthist object containing the telemetry data and then to combine this
with true capture--recapture data (e.g. detections from hair-snag DNA)
in another capthist object.  \code{secr.fit} automatically detects the
telemetry data in the new object.

}

\usage{

addTelemetry (detectionCH, telemetryCH, type = c('concurrent','dependent','independent'), 
   collapsetelemetry = TRUE, verify = TRUE) 

xy2CH (CH, inflation = 1e-08)

telemetrytype (object) <- value

telemetrytype (object, ...)       

}

\arguments{
  \item{detectionCH}{single-session capthist object, detector type `single', `multi', `proximity' or `count'}
  \item{telemetryCH}{single-session capthist object, detector type `telemetryonly'}
  \item{type}{character (see Details)}
  \item{collapsetelemetry}{logical; if TRUE then telemetry occasions are collapsed to one}
  \item{verify}{logical; if TRUE then \code{\link{verify.capthist}} is called on the output}
  \item{CH}{capthist object with telemetryxy attribute}
  \item{inflation}{numeric tolerance for polygon}
  \item{object}{secr traps object}
  \item{value}{character telemetry type replacement value}
  \item{\dots}{other arguments}
}

\details{

It is assumed that a number of animals have been radiotagged, and their telemetry data
(xy-coordinates) have been input to \code{telemetryCH}, perhaps using
\code{\link{read.capthist}} with \code{detector = "telemetryonly"} and \code{fmt =
  "XY"}, or with \code{\link{read.telemetry}}.

A new capthist object is built comprising all the detection
histories in \code{detectionCH}, plus empty (all-zero) histories for
every telemetered animal not in \code{detectionCH}.  Telemetry is associated with new sampling occasions and a new detector (nominally at the same point as the first in \code{detectionCH}).  The number of telemetry fixes of each animal is recorded in the relevant cell of the new capthist object (CH[i, s, K+1] for animal i and occasion s if there were K detectors in detectionCH).

The new sampling occasion(s) are assigned the detector type `telemetry' in the traps attribute of the output capthist object, and the traps attribute \code{telemetrytype} is set to the value provided.  The telemetry type may be ``independent'' (no matching of individuals in captured and telemetered samples), ``dependent'' (telemetered animals are a subset of captured animals) or ``concurrent'' (histories may be capture-only, telemetry-only or both capture and telemetry).

The telemetry locations are carried over from telemetryCH as attribute `xylist' (each
component of xylist holds the coordinates of one animal; use
\code{\link{telemetryxy}} to extract).

The default behaviour of `addTelemetry` is to automatically collapse all telemetry occasions into one.  This is computationally more efficient than the alternative, but closes off some possible models.

\code{xy2CH} partly reverses \code{addTelemetry}: the location
information in the telemetryxy attribute is converted back to a capthist with
detector type `telemetry'.

}

\value{
  A single-session capthist object with the same detector type as
  \code{detectionCH}, but possibly with empty rows and an `telemetryxy' attribute.
}

\note{

Telemetry provides independent data on the location and presence of a
sample of animals. These animals may be missed in the main sampling that
gives rise to detectionCH i.e., they may have all-zero detection
histories.

The `telemetry' detector type is used for telemetry occasions in a combined dataset.

}

\seealso{
  
  \code{\link{capthist}},
  \code{\link{make.telemetry}},
  \code{\link{read.telemetry}},
  \code{\link{telemetryxy}}
  \code{\link{telemetered}}
  
}

\examples{

\dontrun{

# Generate some detection and telemetry data, combine them using
# addTelemetry, and perform analyses

# detectors
te <- make.telemetry()
tr <- make.grid(detector = "proximity")

# simulated population and 50\% telemetry sample
totalpop <- sim.popn(tr, D = 20, buffer = 100)
tepop <- subset(totalpop, runif(nrow(totalpop)) < 0.5)

# simulated detection histories and telemetry
# the original animalID (renumber = FALSE) are needed for matching
trCH <- sim.capthist(tr,  popn = totalpop, renumber = FALSE, detectfn = "HHN")
teCH <- sim.capthist(te, popn = tepop, renumber=FALSE, detectfn = "HHN",
    detectpar = list(lambda0 = 3, sigma = 25))

combinedCH <- addTelemetry(trCH, teCH)

# summarise and display
summary(combinedCH)
plot(combinedCH, border = 150)
ncapt <- apply(combinedCH,1,sum)
points(totalpop[row.names(combinedCH)[ncapt==0],], pch = 1)
points(totalpop[row.names(combinedCH)[ncapt>0],], pch = 16)

# for later comparison of precision we must fix the habitat mask
mask <- make.mask(tr, buffer = 100)
fit.tr <- secr.fit(trCH, mask = mask, CL = TRUE, detectfn = "HHN")  ## trapping alone
fit.te <- secr.fit(teCH, mask = mask, CL = TRUE, start = log(20),   ## telemetry alone
    detectfn = "HHN") 
fit2   <- secr.fit(combinedCH, mask = mask, CL = TRUE,              ## combined
    detectfn = "HHN")                                 

# improved precision when focus on realised population
# (compare CVD)
derived(fit.tr, distribution = "binomial")
derived(fit2, distribution = "binomial")


# may also use CL = FALSE
secr.fit(combinedCH, CL = FALSE, detectfn = "HHN", trace = FALSE)
}

}

\keyword{ manip }
