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

Data from an arbitrary sample of telemetered animals combined with cue counts 
for unmarked animals are a special case that is discussed below 
('Telemetry-calibrated cue counts').

}

\usage{

addTelemetry (detectionCH, telemetryCH, type = c('concurrent','dependent',
   'independent','marking'), collapsetelemetry = TRUE, verify = TRUE, 
   prefix = 'T') 

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
  \item{prefix}{character; pre-pended to ID of telemetry-only animals}
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

The new sampling occasion(s) are assigned the detector type `telemetry' in the 
traps attribute of the output capthist object, and the traps attribute \code{telemetrytype} 
is set to the value provided.  The telemetry type may be ``independent'' 
(no matching of individuals in captured and telemetered samples), ``dependent'' 
(telemetered animals are a subset of captured animals), ``concurrent'' (histories 
may be capture-only, telemetry-only or both capture and telemetry), or ``marking'' 
(one detection history for each telemetered animal, and no others).

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

\section{Telemetry-calibrated cue counts}{

Telemetry type `marking' is used for telemetry-calibrated cue counts: cameras (or other `count' 
detectors) detect \emph{all} animals, but identify only animals wearing 
radio-collars (the `marked' animals). The collared animals also provide 
telemetry fixes. The data are assembled as follows.

\itemize{
\item \code{detectionCH} holds the identified sightings of collared animals on the sighting occasions, one history for each collared animal. Animals never sighted have all-zero histories; they are added by \code{addTelemetry} if they are present in \code{telemetryCH}. Detector type should be `count'.
\item \code{telemetryCH} holds the telemetry fixes. They are collapsed to a single telemetry occasion (the last), recorded at a notional extra detector.
\item The traps attribute \code{markocc} must be set (see \code{\link{markocc}}) with 0 for each sighting occasion and 1 for the telemetry occasion, e.g. \code{c(0,0,0,0,0,1)}.
\item The counts of sightings of animals that were \emph{not} identified as collared are added with \code{\link{addSightings}}: argument \code{unmarked} (attribute Tu, a K x S matrix of counts) for animals without collars and, optionally, \code{nonID} (attribute Tm) for collared animals that were seen but not identified.
}

Cues of collared animals recorded in \code{Tm} or as identified sightings are not counted again as unmarked cues; this is the principal difference from a model that treats the sightings of unmarked animals independently of those of marked animals. The probability that a sighting of a collared animal is identified (parameter \code{pID}, `q') is estimated only if \code{Tm} is present.

A collar may not function throughout the survey. Use \code{marked(capthist) <- } (see \code{\link{marked}}) to specify when each animal could be recognised as marked (default is all occasions). Collared animals that were never detected contribute information only through their fixes and exposure.

The model is fitted by \code{\link{secr.fit}} with \code{detectfn = "HHN"} (or `HEX'). The activity centres of collared animals are integrated by summing over the mask or, usually faster, by Gauss-Hermite quadrature around the mean of the fixes (\code{details = list(telemetryint = "GH")}; see \code{\link{details}}). Density is estimated from the total sightings of unmarked animals after allowing for the collared animals' own sightings, so collars need not be deployed in proportion to density, but the collared animals should be representative of the population in their detectability (absolute density is biased low if collars favour animals that are individually more detectable). A density model such as \code{D ~ x} is possible if the mask has covariate \code{x}; coefficients are then estimated from the spatial pattern of unmarked sightings, and \code{\link{derived}} reports a cue-weighted density `Dcw' (the mean of D(x) weighted by the exposure of the detectors) that is better determined than the mean density over the mask. A mask covariate named `marking' gives a non-uniform prior for the activity centres of collared animals (default uniform).

\code{\link{secr.fit}} stops with an error for telemetry type `marking' if any non-telemetry detector is not of type `count'.

With \code{CL = TRUE} the likelihood is conditional on the collared animals: it uses their identified sightings, their telemetry and, if present, \code{Tm}. The unmarked sightings \code{Tu} contribute only through their spatial pattern across detectors, and only when there is a (relative) density model such as \code{D ~ x}; the density level is not estimated. The detection parameters (including \code{pID}) and the density coefficients are those of the full likelihood, and the CL fit is a quick way to obtain them. Density is obtained by refitting with \code{CL = FALSE}, using the CL fit as \code{start}; \code{\link{derived}}, \code{\link{derivedDcoef}}, \code{\link{derivedDsurface}} and \code{\link{region.N}} are not available for a CL fit of this type.

Sightings are overdispersed relative to the Poisson (the unmarked total is a sum over animals of correlated cues). The scale factor `chat' is computed analytically from the fitted model: fit once, then \code{secr.refit(fit, details = list(chatmethod = "analytic"))}. Simulation of \code{chat} (\code{details$nsim}) is not available for this telemetry type. See \code{\link{details}}.

Finite mixture models (e.g. \code{model = list(sigma ~ h2, lambda0 ~ h2)}) are supported when the class of every collared animal is known: the classes are given by a categorical individual covariate named in \code{hcov}, and an error is reported if any collared animal has an unknown class. The classes of the unmarked animals are unknown and are modelled as a latent mixture with population proportions \code{pmix}. By default \code{pmix} is estimated from the class composition of the collared animals, which assumes that collaring was independent of class. If that is doubtful (for example, sex-biased collaring) fix \code{pmix} at an external value with \code{fixed = list(pmix = p)}, where for two classes \code{p} is the proportion in the second class.

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
  \code{\link{telemetryxy}},
  \code{\link{telemetered}},
  \code{\link{addSightings}},
  \code{\link{marked}},
  \code{\link{markocc}}
  
}

\examples{

# Telemetry type 'marking': simulated cameras count all animals, but only
# collared animals are identified

set.seed(1)
tr  <- make.grid(nx = 9, ny = 9, spacing = 11, detector = "count")
msk <- make.mask(tr, buffer = 30, spacing = 5, type = "traprect")
pop <- sim.popn(D = 50, core = msk, buffer = 0)
collared <- 1:nrow(pop) \%in\% sample(nrow(pop), 15)
dp <- list(lambda0 = 1, sigma = 4)

# identified sightings of collared animals, and telemetry fixes
sightCH <- sim.capthist(tr, popn = subset(pop, collared), detectfn = "HHN",
    detectpar = dp, noccasions = 5, renumber = FALSE)
teCH <- sim.capthist(make.telemetry(), popn = subset(pop, collared),
    detectfn = "HHN", detectpar = dp, noccasions = 50, renumber = FALSE)

# counts at each detector and occasion of animals without collars (K x S matrix)
uCH <- sim.capthist(tr, popn = subset(pop, !collared), detectfn = "HHN",
    detectpar = dp, noccasions = 5, renumber = FALSE)
Tu <- t(apply(uCH, 2:3, sum))

ch <- addTelemetry(sightCH, teCH, type = "marking")
markocc(traps(ch)) <- c(0,0,0,0,0,1)    # 5 sighting occasions + telemetry
ch <- addSightings(ch, unmarked = Tu)

\donttest{
fit <- secr.fit(ch, mask = msk, detectfn = "HHN", trace = FALSE,
    details = list(telemetryint = "GH"))
predict(fit)

# model-based overdispersion (chat) and cue-weighted density
fit2 <- secr.refit(fit, details = list(chatmethod = "analytic"), trace = FALSE)
predict(fit2)
derived(fit2)
}

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

