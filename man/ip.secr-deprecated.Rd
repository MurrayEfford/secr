\name{ip.secr-deprecated}
\alias{ip.secr-deprecated}
\alias{ip.secr}
\alias{pfn}

\title{

Deprecated Functions

}
\description{

Estimate population density by simulation and inverse prediction (Efford 2004; Efford,
Dawson & Robbins 2004).  A restricted range of SECR models may be fitted
(detection functions with more than 2 parameters are not supported, nor
are covariates).

These functions are replaced by \code{ipsecr.fit} and \code{proxyfn1}
in the forthcoming package \pkg{ipsecr} (https://github.com/MurrayEfford/ipsecr). 
They will be removed in future versions of \pkg{secr}.

}

\usage{

ip.secr (capthist, predictorfn = pfn, predictortype = "null", detectfn = 0,
    mask = NULL, start = NULL, boxsize = 0.2, boxsize2 = boxsize, centre = 3,  
    min.nsim = 10, max.nsim = 2000, CVmax = 0.002, var.nsim = 1000, maxbox = 5, 
    maxtries = 2, ncores = 1, seed = NULL, trace = TRUE, ...)

pfn(capthist, N.estimator = c("n", "null","zippin","jackknife") )

}

\arguments{
  \item{capthist}{ capthist object including capture data and detector (trap) layout}
  \item{predictorfn}{ a function with two arguments (the first a capthist object) that returns a vector of predictor values }
  \item{predictortype}{ value (usually character) passed as the second argument of predictorfn }
  \item{detectfn}{integer code or character string for shape of detection function 0 halfnormal, 2 exponential, 3 uniform) -- see \link{detectfn} }
  \item{mask}{ optional habitat mask to limit simulated population }
  \item{start}{vector of np initial parameter values (density, g0 and sigma)}
  \item{boxsize}{ scalar or vector of length np for size of design as fraction of central parameter value}
  \item{boxsize2}{ as for \code{boxsize}; used from second box onwards}
  \item{centre}{ number of centre points in simulation design}
  \item{min.nsim}{ minimum number of simulations per point }
  \item{max.nsim}{ maximum number of simulations per point }
  \item{CVmax}{ tolerance for precision of points in predictor space }
  \item{var.nsim}{ number of additional simulations to estimate variance-covariance matrix }
  \item{maxbox}{ maximum number of attempts to `frame' solution }
  \item{maxtries}{ maximum number of attempts at each simulation }
  \item{ncores}{integer number of cores available for parallel processing}  
  \item{seed}{either NULL or an integer that will be used in a call to \code{set.seed}}
  \item{trace}{logical, set FALSE to suppress progress reports}
  \item{\dots}{ further arguments passed to sim.popn }
  \item{N.estimator}{ character value indicating population estimator to use }
}

\details{

`Inverse prediction' uses methods from multivariate calibration (Brown
1982). The goal is to estimate population density (D) and the parameters
of a detection function (usually g0 and sigma) by `matching' statistics
from \code{predictorfn(capthist)} (the target vector) and statistics
from simulations of a 2-D population using the postulated detection
model. Statistics (see Note) are defined by the predictor function,
which should return a vector equal in length to the number of parameters
(np = 3). Simulations of the 2-D population use \code{\link{sim.popn}}.
The simulated population is sampled with \code{\link{sim.capthist}}
according to the detector type (e.g., `single' or `multi') and detector
layout specified in traps(capthist), including allowance for varying
effort if the layout has a \code{\link{usage}} attribute.

\dots may be used to control aspects of the simulation by passing named
arguments (other than D) to \code{sim.popn}. The most important arguments of
\code{sim.popn} to keep an eye on are `buffer' and `Ndist'. `buffer'
defines the region over which animals are simulated (unless \code{mask}
is specified) - the region should be large enough to encompass all
animals that might be caught. `Ndist' controls the number of individuals
simulated within the buffered or masked area. The default is `poisson'.
Use `Ndist = fixed' to fix the number in the buffered or masked area
\eqn{A} at \eqn{N = DA}. This conditioning reduces the estimated
standard error of \eqn{\hat{D}}{D-hat}, but conditioning is not always
justified - seek advice from a statistician if you are unsure.

The simulated 2-D distribution of animals is Poisson by default. There
is no `even' option as in Density.

Simulations are conducted on a factorial experimental design in
parameter space - i.e. at the vertices of a cuboid `box' centred on the
working values of the parameters, plus an optional number of centre
points. The size of the `box' is specified as a fraction of the working
values, so for example the limits on the density axis are D*(1--boxsize)
and D*(1+boxsize) where D* is the working value of D. For g0, this
computation uses the odds transformation (g0/(1--g0)). \code{boxsize}
may be a vector defining different scaling on each parameter dimension.

A multivariate linear model is fitted to predict each set of simulated
statistics from the known parameter values. The number of simulations at
each design point is increased (doubled) until the residual standard
error divided by the central value is less than CVmax for all
parameters. An error occurs if max.nsim is exceeded.

Once a model with sufficient precision has been obtained, a new working
vector of parameter estimates is `predicted' by inverting the linear
model and applying it to the target vector. A working vector
is accepted as the final estimate when it lies within the box; this
reduces the bias from using a linear approximation to extrapolate a
nonlinear function. If the working vector lies outside the box then a new
design is centred on value for each parameter in the working vector.

Once a final estimate is accepted, further simulations are conducted to
estimate the variance-covariance matrix. These also provide a parametric
bootstrap sample to evaluate possible bias. Set var.nsim = 0
to suppress the variance step.

See Efford et al. (2004) for another description of the method, and
Efford et al. (2005) for an application.

The value of \code{predictortype} is passed as the second argument of
the chosen \code{predictorfn}. By default this is \code{pfn}, for which
the second argument (\code{N.estimator}) is a character value from
c("n", "null","zippin","jackknife"), corresponding respectively to the
number of individuals caught (Mt+1), and \eqn{\hat{N}}{N-hat} from
models M0, Mh and Mb of Otis et al. (1978).

If not provided, the starting values are determined automatically
with \code{\link{autoini}}.

Linear measurements are assumed to be in metres and density in animals
per hectare (10 000 \eqn{\text{m}^2}{m^2}).

If \code{ncores > 1} the \pkg{parallel} package is used to create
processes on multiple cores (see \link{Parallel} for more).

}

\value{
For \code{ip.secr}, a list comprising
  \item{call}{ the function call }
  \item{IP}{ dataframe with estimated density \eqn{\text{ha}^{-1}}{/ha}, g0 and sigma (m) }
  \item{vcov}{ variance-covariance matrix of estimates }
  \item{ip.nsim}{ total number of simulations }
  \item{variance.bootstrap}{ dataframe summarising simulations for variance estimation }
  \item{proctime}{ processor time (seconds) }

For \code{pfn}, a vector of numeric values corresponding to \eqn{\hat{N}}{N-hat}, \eqn{\hat{p}}{p-hat}, and
\code{RPSV}, a measure of the spatial scale of individual detections.

}

\references{

Brown, P. J. (1982) Multivariate calibration. \emph{Journal of the Royal
Statistical Society, Series B} \bold{44}, 287--321.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture--recapture: likelihood-based methods. In:
D. L. Thompson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer. Pp. 255--269.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
\emph{Animal Biodiversity and Conservation} \bold{27},
 217--228.

Efford, M. G., Warburton, B., Coleman, M. C. and Barker, R. J. (2005) A
field test of two methods for density estimation. \emph{Wildlife Society
Bulletin} \bold{33}, 731--738.

Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978)
Statistical inference from capture data on closed animal populations.
\emph{Wildlife Monographs} \bold{62}.

}

\section{Warning}{ Simulation becomes unreliable with very sparse
  populations, or sparse sampling, because some simulated datasets will
  have no recaptures or even no captures. Adjustments were made in secr
  2.3.1 to make the function more stable in these conditions (e.g.,
  allowing a failed simulation to be repeated, by setting the `maxtries'
  argument > 1), but results probably should not be relied upon when
  there are warning messages regarding failed simulations. }

\note{

Each statistic is expected to have a monotonic relationship with one
parameter when the other parameters are held constant. Typical statistics are -

\tabular{ll}{
Statistic \tab Parameter \cr
\eqn{\hat{N}}{N-hat} \tab \eqn{D} \cr
\eqn{\hat{p}}{p-hat} \tab \eqn{g_0}{g0} \cr
\eqn{RPSV} \tab \eqn{\sigma}{sigma}
}

where \eqn{\hat{N}}{N-hat} and \eqn{\hat{p}}{p-hat} are estimates of population size and capture probability from the
naive application of a nonspatial population estimator, and \eqn{RPSV}
is a trap-revealed measure of the scale of movement.

This method provides nearly unbiased estimates of the detection
parameter g0 when data are from single-catch traps (likelihood-based
estimates of g0 are biased in this case - Efford, Borchers & Byrom 2009).

The implementation largely follows that in Density, and it may help to
consult the Density online help. There are some differences: the M0 and
Mb estimates of population-size in \code{ip.secr} can take non-integer
values; the simulation design used by \code{ip.secr} uses odds(g0)
rather than g0; the default boxsize and CVmax differ from those in
Density 4.5. There is no provision in \code{ip.secr} for two-phase
estimation, using a different experimental design at the second phase.
If you wish you can achieve the same effect by using the estimates as
starting values for a second call of \code{ip.secr} (see examples).

Maximum likelihood estimates from \code{secr.fit} are preferable in
several respects to estimates from inverse prediction (speed*; more complex
models; tools for model selection). \code{ip.secr} is provided for checking
estimates of g0 from single-catch traps, and for historical continuity.

* \code{\link{autoini}} with thin = 1 provides fast estimates
from a simple halfnormal model if variances are not required.

}

\seealso{

\code{\link{capthist}}, \code{\link{secr.fit}}, \code{\link{RPSV}}, \code{\link{autoini}},
    \code{\link{sim.popn}}, \code{\link{Detection functions}}

}

\examples{

\dontrun{

## these calculations may take several minutes

## default settings
ip.secr (captdata)

## coarse initial fit, no variance step
ip1 <- ip.secr (captdata, boxsize = 0.2, CVmax=0.01, var=0)
## refined fit
ip2 <- ip.secr (captdata, start = ip1$IP[,"estimate"],
    boxsize = 0.1, CVmax=0.002, var=1000)
ip2

## compare to MLE of same data using multi-catch assumption
predict(secrdemo.0)

## improvise another predictor function (dbar instead of RPSV)
pfn2 <- function (capthist, v) {  ## v is not used
    sumni <- sum(capthist!=0)   ## total detections
    n <- nrow(capthist)         ## number of individuals
    nocc <- ncol(capthist)      ## number of occasions
    c(N = n, p = sumni/n/nocc, dbar = dbar(capthist))
}
ip.secr (captdata, predictorfn = pfn2)
}
}

\keyword{ models }

