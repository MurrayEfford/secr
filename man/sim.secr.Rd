\name{sim.secr}
\alias{simulate.secr}
\alias{sim.secr}
\alias{sim.detect}
\alias{simulate}

\title{
  Simulate From Fitted secr Model
}

\description{

Simulate a spatially distributed population, sample from that population
with an array of detectors, and optionally fit an SECR model to the
simulated data.

}

\usage{

\method{simulate}{secr}(object, nsim = 1, seed = NULL, maxperpoly = 100,
    chat = 1, ...)

sim.secr(object, nsim = 1, extractfn = function(x) c(deviance =
    deviance(x), df = df.residual(x)), seed = NULL, maxperpoly = 100,
    data = NULL, tracelevel = 1, hessian = c("none", "auto", "fdHess"),
    start = object$fit$par, ncores = NULL)

sim.detect(object, popnlist, maxperpoly = 100, renumber = TRUE)

}

\arguments{
  \item{object}{a fitted secr model}
  \item{nsim}{integer number of replicates}
  \item{seed}{either NULL or an integer that will be used in a call to \code{set.seed}}
  \item{maxperpoly}{ integer maximum number of detections of an individual
  in one polygon or transect on any occasion }
  \item{chat}{real value for overdispersion parameter}
  \item{ncores}{integer number of threads used by \code{secr.fit}}
  \item{extractfn}{function to extract output values from fitted model}
  \item{data}{optional list of simulated data saved from previous call to \code{simulate.secr}}
  \item{tracelevel}{integer for level of detail in reporting (0,1,2)}
  \item{hessian}{character or logical controlling the computation of the Hessian matrix}
  \item{start}{vector of starting `beta' values for \code{secr.fit}}
  \item{...}{other arguments (not used)}
  \item{popnlist}{list of popn objects}
  \item{renumber}{logical; if TRUE then }
}

\details{

For each replicate, \code{simulate.secr} calls \code{\link{sim.popn}} to
generate session- and group-specific realizations of the (possibly
inhomogeneous) 2-D Poisson distribution fitted in \code{object}, across
the habitat mask(s) in \code{object}.  Group subpopulations are combined
using \code{\link{rbind.popn}} within each session; information to
reconstruct groups is retained in the individual-level factor
covariate(s) of the resulting \code{popn} object (corresponding to
\code{object$groups}).  Each population is then sampled using the fitted
detection model and detector (trap) array(s) in \code{object}.

The random number seed is managed as in \code{simulate.lm}.

Certain model types are not supported by \code{simulate.secr}. These
include models fitted using conditional likelihood (\code{object$CL =
TRUE}), telemetry models and exotic behavioural response models.

Detector type is determined by \code{detector(traps(object$capthist))}.

\code{sim.secr} is a wrapper function. If \code{data = NULL} (the
default) then it calls \code{simulate.secr} to generate \code{nsim} new datasets. If
\code{data} is provided then \code{nsim} is taken to be
\code{length(data)}. \code{secr.fit} is called to fit the original model
to each new dataset. Results are summarized according to the
user-provided function \code{extractfn}. The default \code{extractfn}
returns the deviance and its degrees of freedom; a NULL value for
\code{extractfn} returns the fitted secr objects after
\code{\link{trim}}ming to reduce bulk. Simulation uses the detector type
of the data, even when another likelihood is fitted (this is the case
with single-catch data, for which a multi-catch likelihood is fitted).
Warning messages from \code{secr.fit} are suppressed.

\code{extractfn} should be a function that takes an \code{secr} object
as its only argument.

\code{tracelevel=0} suppresses most messages; \code{tracelevel=1} gives a
terse message at the start of each fit; \code{tracelevel=2} also sets
`details$trace = TRUE' for \code{secr.fit}, causing each likelihood
evaluation to be reported.

\code{hessian} controls computation of the Hessian matrix from which
variances and covariances are obtained. \code{hessian} replaces the
value in \code{object\$details}.  Options are "none" (no variances),
"auto" (the default) or "fdhess" (see \code{\link{secr.fit}}). It is OK
(and faster) to use \code{hessian="none"} unless \code{extractfn} needs
variances or covariances. Logical TRUE and FALSE are interpreted by
\code{secr.fit} as "auto" and "none".

If \code{ncores = NULL} then the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS is used (see \code{\link{setNumThreads}}).

\code{sim.capthist} is a more direct way to simulate data from a null
model (i.e. one with constant parameters for density and detection), or
from a time-varying model.

\code{sim.detect} is a function used internally that will not usually be
called directly.

}
\value{

For \code{simulate.secr}, a list of data sets (`capthist' objects). This
list has class \code{c("secrdata", "list")}; the initial state of the random
number generator (roughly, the value of .Random.seed) is stored as the
attribute `seed'.

The value from \code{sim.secr} depends on \code{extractfn}: if that
returns a numeric vector of length \code{n.extract} then the value is a
matrix with \code{dim = c(nsim, n.extract)} (i.e., the matrix has one
row per replicate and one column for each extracted value). Otherwise,
the value returned by \code{sim.secr} is a list with one component per
replicate (strictly, an object of class = \code{c("secrlist", "list")}). Each
simulated fit may be retrieved \emph{in toto} by specifying
\code{extractfn = identity}, or slimmed down by specifying
\code{extractfn = NULL} or \code{extractfn = trim}, which are
equivalent.

For either form of output from \code{sim.secr} the initial state of the
random number generator is stored as the attribute `seed'.

For \code{sim.detect} a list of `capthist' objects.

}

\seealso{
\code{\link{sim.capthist}}, \code{\link{secr.fit}},
\code{\link[stats]{simulate}}, \code{\link{secr.test}}
}

\note{

The value returned by \code{simulate.secr} is a list of `capthist'
objects; if there is more than one session, each `capthist' is itself a
sort of list .

The classes `secrdata' and `secrlist' are used only to override the ugly
and usually unwanted printing of the seed attribute. However, a few
other methods are available for `secrlist' objects
(e.g. \code{plot.secrlist}).

The default value for \code{start} in \code{sim.secr} is the previously
fitted parameter vector. Alternatives are NULL or \code{object$start}.

}

\section{Warning}{

\code{sim.secr} does not work for mark--resight models.

}

\examples{

\dontrun{

## previously fitted model
simulate(secrdemo.0, nsim = 2)

## The following has been superceded by secr.test()

## this would take a long time...
sims <- sim.secr(secrdemo.0, nsim = 99)
deviance(secrdemo.0)
devs <- c(deviance(secrdemo.0),sims$deviance)
quantile(devs, probs=c(0.95))
rank(devs)[1] / length(devs)

## to assess bias and CI coverage
extrfn <- function (object) unlist(predict(object)["D",-1])
sims <- sim.secr(secrdemo.0, nsim = 50, hessian = "auto",
    extractfn = extrfn)
sims

## with a larger sample, could get parametric bootstrap CI
quantile(sims[,1], c(0.025, 0.975))

}

}
\keyword{ datagen }
\keyword{ models }

