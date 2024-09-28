\name{ovensong}
\docType{data}

\alias{ovensong}
\alias{signalCH}
\alias{ovensong.model.1}
\alias{ovensong.model.2}

\title{ Ovenbird Acoustic Dataset }

\description{

Data from an acoustic survey of ovenbirds (\emph{Seiurus
aurocapilla}) at a site in Maryland, USA.

}

\usage{

signalCH
ovensong.model.1
ovensong.model.2

}

\details{

In June 2007 D. K. Dawson and M. G. Efford used a moving 4-microphone array
to survey breeding birds in deciduous forest at the Patuxent Research
Refuge near Laurel, Maryland, USA. The data for ovenbirds were used to
demonstrate a new method for analysing acoustic data (Dawson and Efford
2009). See \code{\link{ovenbird}} for mist-netting data from the same
site over 2005--2009, and for other background.

Over five days, four microphones were placed in a square (21-m side)
centred at each of 75 points in a rectangular grid (spacing 50 m); on
each day points 100 m apart were sampled sequentially. Recordings of 5
minutes duration were made in .wav format on a 4-channel digital sound
recorder.

The data are estimates of average power on each channel (microphone) for
the first song of each ovenbird distinguishable in a particular 5-minute
recording. Power was estimated with the sound analysis software Raven
Pro 1.4 (Charif et al. 2008), using a window of 0.7 s duration and
frequencies between 4200 and 5200 Hz, placed manually at the approximate
centre of each ovenbird song. Sometimes this frequency range was
obscured by insect noise so an alternative 1000-Hz range was measured
and the values were adjusted by regression.

The data are provided as a single-session, single-occasion
\code{capthist} object \code{signalCH}. The `signal' attribute contains
the power measurement in decibels for each detected sound on each
channel where the power threshold is exceeded. As the threshold signal
(attribute \code{cutval} = 35) is less than any signal value in this
dataset, all detection histories are complete (1,1,1,1) across
microphones. For analysis Dawson and Efford applied a higher threshold
that treated weaker signals as `not detected' (see Examples).

The row names of \code{signalCH} (e.g. "3755AX") are formed from a
4-digit number indicating the sampling location (one of 75 points on a
50-m grid) and a letter A--D to distinguish individual ovenbirds within
a 5-minute recording; `X' indicates power values adjusted by regression.

The default model for sound attenuation is a log-linear decline with
distance from the source (linear decline on dB scale). Including a
spherical spreading term in the sound attenuation model causes the
likelihood surface to become multimodal in this case. Newton-Raphson,
the default maximization method in \code{secr.fit}, is particularly
inclined to settle on a local maximum; in the example below we use a set
of starting values that have been found by trial and error to yield the
global maximum.

Two fitted models are included (see Examples for details).

\tabular{ll}{
  Object \tab Description \cr
  signalCH \tab capthist object \cr
  ovensong.model.1 \tab fitted secr model -- spherical spreading \cr
  ovensong.model.2 \tab fitted secr model -- no spherical spreading \cr
}

}

\source{

D. K. Dawson (\email{ddawson@usgs.gov}) and M. G. Efford unpublished data.

}

\references{

Charif, R. A., Waack, A. M. and Strickman, L. M. (2008) Raven Pro 1.3 User's
Manual. Cornell Laboratory of Ornithology, Ithaca, New York.

Dawson, D. K. and Efford, M. G. (2009) Bird population density estimated from
acoustic signals. \emph{Journal of Applied Ecology} \bold{46}, 1201--1209.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population density estimated
from locations of individuals on a passive detector array. \emph{Ecology} \bold{90},
2676--2682.

}

\seealso{

\code{\link{capthist}}, \code{\link{ovenbird}}, \code{\link{Detection functions}} 

}

\examples{

summary(signalCH)
traps(signalCH)
signal(signalCH)

## apply signal threshold
signalCH.525 <- subset(signalCH, cutval = 52.5)

\dontrun{

## models with and without spherical spreading
omask <- make.mask(traps(signalCH), buffer = 200)
ostart <- c(log(20), 80, log(0.1), log(2))
ovensong.model.1 <- secr.fit( signalCH.525, mask = omask, 
    start = ostart, detectfn = 11 ) 
ovensong.model.2 <- secr.fit( signalCH.525, mask = omask, 
    start = ostart, detectfn = 10 ) 

}

## compare fit of models
AIC(ovensong.model.1, ovensong.model.2)

## density estimates, dividing by 75 to allow for replication
collate(ovensong.model.1, ovensong.model.2)[1,,,"D"]/75

## plot attenuation curves cf Dawson & Efford (2009) Fig 5
pars1 <- predict(ovensong.model.1)[c("beta0", "beta1"), "estimate"]
pars2 <- predict(ovensong.model.2)[c("beta0", "beta1"), "estimate"]
attenuationplot(pars1, xval=0:150, spherical = TRUE, ylim = c(40,110))
attenuationplot(pars2, xval=0:150, spherical = FALSE, add = TRUE, 
    col = "red")
## spherical spreading only
pars1[2] <- 0  
attenuationplot(pars1, xval=0:150, spherical = TRUE, add = TRUE, lty=2)

}

