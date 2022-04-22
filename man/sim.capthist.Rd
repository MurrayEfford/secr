\name{sim.capthist}
\alias{sim.capthist}
\alias{sim.resight}

\title{ Simulate Detection Histories }

\description{

Create a set of capture or marking-and-resighting histories by simulated sampling
of a 2-D population using an array of detectors.

}

\usage{
sim.capthist(traps, popn = list(D = 5, buffer = 100,
    Ndist = "poisson"), detectfn = 0, detectpar = list(),
    noccasions = 5, nsessions = 1, binomN = NULL, exactN = NULL,
    p.available = 1, renumber = TRUE, seed = NULL,
    maxperpoly = 100, chulltol = 0.001, userdist = NULL, 
    savepopn = FALSE)
sim.resight(traps, popn = list(D = 5, buffer = 100, Ndist = "poisson"), ...,
    pID = 1, unmarked = TRUE, nonID = TRUE,  unresolved = FALSE, unsighted = TRUE, 
    pmark = 0.5, Nmark = NULL, markingmask = NULL)
}

\arguments{
  \item{traps}{ \code{traps} object with the locations and other attributes of detectors }
  \item{popn}{ locations of individuals in the population to be sampled, either as
 a \code{\link{popn}} object (see \code{\link{sim.popn}}) or a list with the named components. }
  \item{detectfn}{integer code or character string for shape of
  detection function 0 = halfnormal etc. -- see \link{detectfn}}
  \item{detectpar}{ list of values for named parameters of detection function }
  \item{noccasions}{ number of occasions to simulate }
  \item{nsessions}{ number of sessions to simulate }
  \item{binomN}{ integer code for distribution of counts (see Details) }
  \item{exactN}{ integer number of telemetry fixes per occasion }
  \item{p.available}{ vector of one or two probabilities (see Details) }
  \item{renumber}{ logical for whether output rows should labeled sequentially (TRUE) or
 retain the numbering of the population from which they were drawn (FALSE) }
  \item{seed}{either NULL or an integer that will be used in a call to \code{set.seed}}
  \item{maxperpoly}{ integer maximum number of detections of an individual
    in one polygon or transect on any occasion }
  \item{chulltol}{numeric buffer (m) for polygon around telemetry
    locations} 
  \item{userdist}{user-defined distance function or matrix (see \link{details})}  
  \item{savepopn}{ logical; if TRUE then the popn (input or simulated) is saved as an attribute}
  \item{\dots}{ arguments to pass to \code{sim.capthist} }
  \item{pID}{ probability of individual identification for marked animals }
  \item{unmarked}{ logical, if TRUE unmarked individuals are not recorded during `sighting' }
  \item{nonID}{logical, if TRUE then unidentified marked individuals are not recorded during `sighting' }
  \item{unresolved}{logical, if TRUE then individuals of unresolved mark status are not recorded during `sighting' }
  \item{unsighted}{logical, if TRUE and sighting only then capthist includes all-zero histories }
  \item{pmark}{ numeric probability that an individual is `pre-marked' (see Details)}
  \item{Nmark}{ number of individuals to be `pre-marked' (see Details)}
  \item{markingmask}{ \code{mask} object}

}
\details{

If \code{popn} is not of class `popn' then a homogeneous Poisson
population with the desired density (animals/ha) is first simulated over
the rectangular area of the bounding box of \code{traps} plus a buffer
of the requested width (metres). The detection algorithm depends on the
detector type of \code{traps}. For `proximity' detectors, the actual
detection probability of animal \emph{i} at detector \emph{j} is the
naive probability given by the detection function. For `single' and
`multi' detectors the naive probability is modified by competition
between detectors and, in the case of `single' detectors, between animals. See
Efford (2004) and other papers below for details.

Detection parameters in \code{detectpar} are specific to the detection
function, which is indicated by a numeric code (\code{\link{detectfn}}).
Parameters may vary with time - for this provide a vector of length
\code{noccasions}.  The g0 parameter may vary both by time and detector
- for this provide a matrix with \code{noccasions} rows and as many
columns as there are detectors. The default detection parameters are
\code{list(g0 = 0.2, sigma = 25, z = 1)}.

The default is to simulate a single session. This may be overridden by
providing a list of populations to sample (argument \code{popn}) or by
specifying \code{nsessions} > 1 (if both then the number of sessions must
match). Using \code{nsessions} > 1 results in replicate samples of
populations with the same density etc. as specified directly in the
\code{popn} argument.

\code{binomN} determines the statistical distribution of the number of
detections of an individual at a particular `count' detector or polygon
on a particular occasion. A Poisson distribution is indicated by
\code{binomN = 0}; see \code{\link{secr.fit}} for more. The distribution
is always Bernoulli (binary) for `proximity' and `signal' detectors.

If \code{exactN} is not specified or zero then the number of telemetry
fixes is a random variable determined by the other detection settings.

\code{p.available} specifies temporary non-availability for detection in
multi-session simulations. If a single probability is specified then
temporary non-availability is random (independent from session to
session). If two probabilities are given then non-availability is
Markovian (dependent on previous state) and the two values are for
animals available and not available at the preceding session. In the
Markovian case, availability in the first session is assigned at random
according to the equilibrium probability p2 / (1 - p1 + p2). Incomplete
availability is not implemented for sampling lists of populations.

\code{detectpar} may include a component `truncate' for the distance
beyond which detection probability is set to zero. By default this value
is NULL (no specific limit).

\code{detectpar} may also include a component `recapfactor' for a
general learned trap response. For `single' and `multi' detector types
the probability of detection changes by this factor for all occasions
after the occasion of first capture. Attempted use with other detector
types causes an error. If recapfactor x g(d) > 1.0, g(d) is truncated at
1.0. Other types of response (site-specific bk, Markovian B) are not
allowed.

If \code{popn} is specified by an object of class `popn' then any
individual covariates will be passed on; the \code{covariates} attribute
of the output is otherwise set to NULL.

The random number seed is managed as in \code{\link[stats]{simulate}}.

\code{chulltol} is used only when simulating telemetry locations. By
default, a new 'traps' polygon is generated as the convex hull of the
simulated locations, with a slight (1 mm) added buffer to ensure
boundary points are within the polygon. Buffering is suppressed if
\code{chulltol} is NA or negative.

\code{userdist} cannot be set if `traps' is any of polygon, polygonX,
transect or transectX.

\code{sim.resight} generates mark-resight data. The `markocc' attribute 
of `traps' indicates the occasions which are for sighting-only (0) or 
marking and recapture (1). The number 
of occasions is determined by \code{markocc}. \code{sim.capthist} is first 
called with the arguments `traps' and \dots. The detector type of `traps' 
should be `proximity' or `count' for sighting occasions (markocc = 0). The detector type need not be the same for marking and sighting occasions ('multi' is allowed on marking occasions). If \dots includes a non-null `seed' the 
random seed is reset in \code{sim.resight} and not passed to 
\code{sim.capthist}.

A special case arises when all occasions are sighting-only. Then it is assumed that individuals in the population are marked prior to the start of sampling with a known spatial distribution (i.e. marking does not follow a spatial detection model). By default, animals throughout the buffered area are pre-marked with probability \code{pmark}. If \code{Nmark} is specified then a sample of size \code{Nmark} will be selected for marking, overriding \code{pmark}.

The marked population may be restricted to a subset of the space spanned by \code{popn} by specifying \code{markingmask}, which may have a further covariate `marking' to vary the intensity of marking.

}
\value{

For \code{sim.capthist}, an object of class \code{capthist}, a 3-dimensional array 
with additional attributes. Rows represent
individuals and columns represent occasions; the third dimension, codes the number of 
detections at each detector (zero or one for trap detectors (`single', `multi') and binary proximity detectors.

The initial state of the R random number generator is stored in the
`seed' attribute.

For \code{sim.resight}, an object of class \code{capthist} for which the traps object has a markocc attribute (marking occasions), and there are further attributes Tu (sightings of unmarked animals) and Tm (sightings of marked but not 
identified animals).

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\note{

External code is called to speed the simulations. The present version
assumes a null model, i.e., naive detection probability is constant
except for effects of distance and possibly time (using vector-valued
detection parameters from 1.2.10). You can, however, use
\code{\link{rbind.capthist}} to combine detections of population
subclasses (e.g. males and females) simulated with different parameter
values. This is not valid for detector type "single" because it fails to
allow for competition for traps between subclasses. Future versions may
allow more complex models.

\code{truncate} has no effect (i) when using a uniform detection
function with radius (\code{sigma}) <= \code{truncate} and (ii) with
signal strength detection (detectfn 10, 11). Note that truncated
detection functions are provided for de novo simulation, but are not
available when fitting models with in \code{secr.fit} or simulating from
a fitted model with \code{sim.secr}.

\code{maxperpoly} limits the size of the array allocated for
detections in C code; an error results if the is number is exceeded.

Prior to 2.10.0 \code{sim.resight} interpreted length-2 vectors of detection parameters as referring to marking and sighting occasions; this feature has been discontinued.

}

\seealso{
    \code{\link{sim.popn}}, \code{\link{capthist}}, \code{\link{traps}},
    \code{\link{popn}}, \code{\link{Detection functions}}, \code{\link[stats]{simulate}},
    \code{\link{Tu}}, \code{\link{Tm}}
}

\examples{
## simple example
## detector = "multi" (default)
temptrap <- make.grid(nx = 6, ny = 6, spacing = 20)
sim.capthist (temptrap, detectpar = list(g0 = 0.2, sigma = 20))

## with detector = "proximity", there may be more than one
## detection per individual per occasion
temptrap <- make.grid(nx = 6, ny = 6, spacing = 20, detector =
    "proximity")
summary(sim.capthist (temptrap, detectpar = list(g0 = 0.2, 
    sigma = 20)))

## marking on occasions 1, 3 only
temptrap <- make.grid(nx = 6, ny = 6, spacing = 20, detector = 'proximity')
markocc(temptrap) <- c(1,0,1,0,0)
CH <- sim.resight (temptrap, detectpar = list(g0 = 0.2, sigma = 20))
summary(CH)

## multiple sessions
grid4 <- make.grid(nx = 2, ny = 2)
temp <- sim.capthist (grid4, popn = list(D = 1), nsessions = 20)
summary(temp, terse = TRUE)

## unmarked or presence types
# grid <- make.grid(nx = 10, ny = 10, detector = "unmarked")
# CH <- sim.capthist (grid, noccasions = 5)
# CH
## "presence" and "unmarked" data are stored as "count" data
## behaviour is controlled by detector type, e.g.
# detector(traps(CH)) <- "presence"
# CH


}
\keyword{ datagen }

