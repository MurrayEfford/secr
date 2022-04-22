\name{expected.n}
\alias{expected.n}

\title{
Expected Number of Individuals
}

\description{

  Computes the expected number of individuals detected across a detector
  layout or at each cluster of detectors. 

}

\usage{

expected.n(object, session = NULL, group = NULL, bycluster
    = FALSE, splitmask = FALSE, ncores = NULL)

}

\arguments{
  \item{object}{ \code{secr} object output from \code{secr.fit} }
  \item{session}{ character session vector }
  \item{group}{ group -- for future use }
  \item{bycluster}{ logical to output the expected number for
  clusters of detectors rather than whole array }
  \item{splitmask}{logical for computation method (see Details)}
  \item{ncores}{ integer number of threads to be used for parallel processing}
}

\details{

  The expected number of individuals detected is \eqn{E(n) = \int p.(X)
  D(X) dX}{E(n) = integral p.(X) D(X) dX} where the integration is a
  summation over \code{object$mask}. \eqn{p.(X)} is the probability an
  individual at \eqn{X} will be detected at least once either on the
  whole detector layout (\code{bycluster = FALSE}) or on the detectors
  in a single cluster (see \link{pdot} for more on \eqn{p.}). \eqn{D(X)}
  is the expected density at \eqn{X}, given the model. \eqn{D(X)} is
  constant (i.e. density surface flat) if \code{object$CL == TRUE} or
  \code{object$model$D == ~1}, and for some other possible models.

If the \code{bycluster} option is selected and detectors are not, in
fact, assigned to clusters then each detector will be treated as a
cluster, with a warning.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

By default, a full habitat mask is used for each cluster. This is the
more robust option. Alternatively, the mask may be split into subregions
defined by the cells closest to each cluster.

The calculation takes account of any fitted continuous model for spatial
variation in density (note Warning).

}

\section{Warning}{

  This function changed slightly between 2.1.0 and 2.1.1, and now
  performs as indicated here when bycluster = TRUE and
  clusters are not specified.

  Clusters of detectors are assumed to be independent (always true with detector types
  `proximity', `count' etc.). The computed E(n) does not apply when
  there is competition among clusters of detectors.
  
The prediction of density at present considers only the base level of
density covariates, such as cell-specific habitat variables.

}

\value{

The expected count (bycluster = FALSE) or a vector of expected counts,
one per cluster. For multi-session data, a list of such vectors.

}

\seealso{
\code{\link{region.N}}
}

\examples{

\dontrun{

expected.n(secrdemo.0)
expected.n(secrdemo.0, bycluster = TRUE)
expected.n(ovenbird.model.D)

## Clustered design
mini <- make.grid(nx = 3, ny = 3, spacing = 50, detector =
    "proximity")
tempgrids <- trap.builder (cluster = mini , method = "all",
    frame = expand.grid(x = seq(1000, 9000, 2000),
    y = seq(1000, 9000, 2000)), plt = TRUE)
capt <- sim.capthist(tempgrids, popn = list(D = 2))
tempmask <- make.mask(tempgrids, buffer = 100,
    type = "clusterbuffer")
fit <- secr.fit(capt, mask = tempmask, trace = FALSE)
En <- expected.n(fit, bycluster = TRUE)

## GoF or overdispersion statistic
p <- length(fit$fit$par)
y <- cluster.counts(capt)
## scaled by n-p
sum((y - En)^2 / En) / (length(En)-p)
sum((y - En)^2 / En) / sum(y/En)

}

}

\keyword{ models }
