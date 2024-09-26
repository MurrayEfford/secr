\name{secr.test}
\alias{secr.test}

\title{Goodness-of-Fit Test}

\description{

  Simple Monte-Carlo goodness-of-fit tests for full-likelihood SECR
  models. The approach is to calculate a statistic from either the raw
  data or a fitted model, and to relate this to the distribution of the
  statistic under the original fitted model. The distribution is
  estimated by simulating data from the model, and possibly re-fitting
  the model to each simulated dataset.

  The suitability of different test statistics has yet to be assessed. 
  An alternative and possibly more informative set of tests is provided 
  in function \code{\link{MCgof}} that implements the approach of 
  Choo et al. (2024) for a range of models.
 
}

\usage{

secr.test(object, nsim = 99, statfn, fit = FALSE, seed =
NULL, ncores = NULL, tracelevel = 1, ...)

}

\arguments{
  
  \item{object}{a fitted secr model}
  
  \item{nsim}{integer number of replicates}

  \item{statfn}{function to compute a numeric vector of one or more statistics from a
    single-session `capthist' object or from a fitted model (see Details)}

  \item{fit}{logical; if TRUE the model is re-fitted to each simulated
  dataset}
 
  \item{seed}{either NULL or an integer that will be used in a call to \code{set.seed}}
  
  \item{ncores}{integer number of threads for parallel processing}
  
  \item{tracelevel}{ see \code{\link{sim.secr}} }
  
  \item{\dots}{other arguments passed to statfn, if needed}
  
}

\details{

  The test statistic(s) may be computed either on a dataset or on a
  fitted model, as determined by the argument \code{fit}. The single
  argument expected by \code{statfn} should be either a `capthist' object
  (\code{fit = FALSE}) or an `secr' object (\code{fit = TRUE}).
  
The default statistic when \code{fit = FALSE} is the proportion of
individuals observed on only one occasion, which is equivalent to
\code{statfn = function(CH) c(f1 = sum(apply(abs(CH) > 0,1,sum) == 1) /
nrow(CH))}. Repeat detections on one occasion at the same or different
detectors are not counted. The default statistic is therefore not
appropriate for some data, specifically from `count' or `polygon'
detectors with few occasions or only one.

The default statistic when \code{fit = TRUE} is the deviance divided by
the residual degrees of freedom (i.e., \code{statfn = function(object)
  c(devdf = deviance(object) / df.residual(object))}). 

The reported probability (p) is the rank of the observed value in the
vector combining the observed value and simulated values, divided by
(nsim + 1). Ranks are computed with \code{\link{rank}} using the default
\code{ties.method = "average"}.

Simulations take account of the usage attribute of detectors in the
original capthist object, given that usage was defined and ignoreusage was not
set.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

\code{statfn} may return a vector of statistics for each observed or
simulated dataset or model: then the value of each statistic will be
calculated for every simulated dataset, and summarised. If \code{fit =
TRUE} the vector of statistics may include both those computed on the
raw data (object$capthist) and the fitted model. 

}

\value{
    An object of class \code{c('secrtest','list')} with components:

    \item{object}{as input}
    \item{nsim}{as input}
    \item{statfn}{as input}
    \item{fit}{as input}
    \item{seed}{as input}
    \item{output}{list comprising the simulated values,
      the observed value(s), and estimated probabilities}

    For multi-session input when \code{fit = FALSE}, `output' is a list
    in which each session provides one component.

    Print and plot methods are provided for `secrtest' objects.
  }

\references{

Choo, Y. R., Sutherland, C. and Johnston, A. (2024) A Monte Carlo resampling 
framework for implementing goodness-of-fit tests in spatial capture-recapture model
\emph{Methods in Ecology and Evolution} DOI: 10.1111/2041-210X.14386.

}

\seealso{
    
  \code{\link{MCgof}},
  \code{\link{print.secrtest}},
  \code{\link{plot.secrtest}},
  \code{\link{simulate.secr}},
  \code{\link{sim.secr}},
  \code{\link{deviance.secr}}
  
}

\note{

  \code{simulate.secr} is always used to simulate the raw data, but
  \code{simulate.secr} does not work for all types of fitted
  model. Models fitted by maximizing the likelihood conditional on
  \eqn{n} (\code{CL = TRUE} in \code{secr.fit}) potentially include
  individual covariates whose distribution in the population is
  unknown. This precludes simulation, and conditional-likelihood models
  in general are therefore not covered by \code{secr.test}.

  Other exclusions include exotic non-binary behavioural responses
  (``bn", ``bkn", ``bkc", ``Bkc" - but these are generally undocumented
  in any case).

  If \code{fit = TRUE} then \code{sim.secr} is used.

  At each simulation a new population is generated across the extent of
  the original mask. If the extent is unduly large then time will be
  wasted simulating the possibility of detection for many essentially
  undetectable animals. This is an argument for keeping the mask tight -
  large enough only to avoid mask-induced bias.
  
}

\examples{

\dontrun{

secr.test(secrdemo.0, nsim = 99)

secr.test(ovenbird.model.1, nsim = 20)

## example combining raw data summary and model fit
## assumes single-session
bothfn <- function(object) {
    CH <- object$capthist
    f1 <- sum(apply(abs(CH) > 0, 1, sum) == 1) / nrow(CH)
    devdf <- deviance(object) / df.residual(object)
    c(f1 = f1, devdf = devdf)
}
test <- secr.test (secrdemo.0, nsim = 19, statfn = bothfn, fit = TRUE)
test
plot(test, main = '')

}

}

\keyword{ htest }
