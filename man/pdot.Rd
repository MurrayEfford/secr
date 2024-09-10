\name{pdot}
\alias{pdot}
\alias{CVpdot}

\title{ Net Detection Probability }

\description{

Compute spatially explicit net probability of detection for individual(s) 
at given coordinates (pdot).

}

\usage{
pdot(X, traps, detectfn = 0, detectpar = list(g0 = 0.2,
    sigma = 25, z = 1), noccasions = NULL, binomN = NULL,
    userdist = NULL, ncores = NULL) 

CVpdot(..., conditional = FALSE)    

}

\arguments{
  \item{X}{ vector or 2-column matrix of coordinates }
  \item{traps}{ \code{traps} object }
  \item{detectfn}{ integer code for detection function q.v. }
  \item{detectpar}{ a named list giving a value for each parameter of detection function }
  \item{noccasions}{ number of sampling intervals (occasions) }
  \item{binomN}{ integer code for discrete distribution (see
    \code{\link{secr.fit}}) }
  \item{userdist}{user-defined distance function or matrix (see \link{userdist})}
  \item{ncores}{integer number of threads}
  
  \item{\dots}{arguments passed to \code{pdot}}
  \item{conditional}{logical; if TRUE then computed mean and CV are conditional on detection}
  
}
\details{

  If \code{traps} has a \link{usage} attribute then \code{noccasions} is
  set accordingly; otherwise it must be provided.
  
The probability computed is \eqn{p.(\mathbf{X}) = 1 - \prod\limits _{k}
\{1 - p_s(\mathbf{X},k)\}^{S}}{p.(X) = 1 - (1 - prod(p_s(X,k))^S} where
the product is over the detectors in \code{traps}, excluding any not
used on a particular occasion. The per-occasion detection function
\eqn{p_s}{p_s} is halfnormal (0) by default, and is assumed not to vary
over the \eqn{S} occasions.

From 4.6.11, the detection parameters g0, lambda0 and sigma for point detectors 
may be detector- and occasion-specific. This is achieved by providing a vector 
of values that is replicated internally to fill a matrix with dimensions 
ntraps x noccasions (i.e. in trap order for occasion 1, then occasion 2 etc.)

For detection functions (10) and (11) the signal threshold `cutval' should be
included in \code{detectpar}, e.g., \code{detectpar = list(beta0 = 103, beta1
= -0.11, sdS = 2, cutval = 52.5)}.

The calculation is not valid for single-catch traps because
\eqn{p.(\mathbf{X})}{p.(X)} is reduced by competition between animals.

\code{userdist} cannot be set if `traps' is any of polygon, polygonX,
transect or transectX. if \code{userdist} is a function requiring
covariates or values of parameters `D' or `noneuc' then \code{X} must
have a covariates attribute with the required columns.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

\code{CVpdot} returns the expected mean and CV of pdot across the points listed in \code{X}, assuming uniform population density. \code{X} is usually a habitat mask. See Notes for details.

}

\value{

For \code{pdot}, a vector of probabilities, one for each row in X.

For \code{CVpdot}, a named vector with elements `meanpdot' and `CVpdot'.

}

\seealso{ 
  \code{\link{secr}}, 
  \code{\link{make.mask}},
  \code{\link{Detection functions}}, 
  \code{\link{pdot.contour}},
  \code{\link{CV}}}

\note{

\code{CVpdot} computes the mean \eqn{\mu} and variance \eqn{V} of the location-specific overall detection probability \eqn{p.(\mathbf{X})}{p.(X)} as follows.

\deqn{\mu = \int p.(\mathbf{X}) f(\mathbf{X}) d\mathbf{X},}{\mu = \int p.(X) f(X) dX, }
\deqn{V = \int p.(\mathbf{X})^2 f(\mathbf{X}) d\mathbf{X} - \mu^2.}{V = \int p.(X)^2 f(X) dX - \mu^2. }

For uniform density and \code{conditional = FALSE}, \eqn{f(\mathbf{X})}{f(X)} is merely a scaling factor independent of \eqn{\mathbf{X}}{X}. 

If \code{conditional = TRUE} then \eqn{f(\mathbf{X}) = p.(\mathbf{X}) / \int p.(\mathbf{X}) d\mathbf{X}}{f(X) = p.(X) /  \int p.(X) dX}. 

The coefficient of variation is CV = \eqn{\sqrt{V}/\mu}{sqrt(V)/\mu}.
}

\examples{

\dontrun{

  temptrap <- make.grid()
  ## per-session detection probability for an individual centred
  ## at a corner trap. By default, noccasions = 5.
  pdot (c(0,0), temptrap, detectpar = list(g0 = 0.2, sigma = 25),
    noccasions = 5)
    
  msk <- make.mask(temptrap, buffer = 100)
  CVpdot(msk, temptrap, detectpar = list(g0 = 0.2, sigma = 25),
    noccasions = 5)
}

}

\keyword{ manip }

