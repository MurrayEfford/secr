\name{region.N}
\alias{region.N}
\alias{region.N.secr}
\alias{region.N.secrlist}
\alias{population size}

\title{Population Size}

\description{

Estimate the expected and realised populations in a region,
using a fitted spatially explicit capture--recapture model. Density is
assumed to follow an inhomogeneous Poisson process in two
dimensions. Expected \eqn{N} is the volume under a fitted density
surface; realised \eqn{N} is the number of individuals within the region
for the current realisation of the process (cf Johnson et al. 2010; see
Note). 

}


\usage{

region.N(object, ...)

\method{region.N}{secr}(object, region = NULL, spacing = NULL, session = NULL,
    group = NULL, se.N = TRUE, alpha = 0.05, loginterval = TRUE,
    keep.region = FALSE, nlowerbound = TRUE, RN.method = "poisson",
    pooled.RN = FALSE, ncores = NULL, ...)
    
\method{region.N}{secrlist}(object, region = NULL, spacing = NULL, session = NULL,
    group = NULL, se.N = TRUE, alpha = 0.05, loginterval = TRUE,
    keep.region = FALSE, nlowerbound = TRUE, RN.method = "poisson",
    pooled.RN = FALSE, ncores = NULL, ...)

}

\arguments{
  
  \item{object}{ \code{secr} object output from \code{secr.fit} }
  
  \item{region}{ mask object defining the possibly non-contiguous region
  for which population size is required, or vector polygon(s) (see
  Details)}

  \item{spacing}{ spacing between grid points (metres) if region mask is
  constructed on the fly }

  \item{session}{ character session }
  
  \item{group}{ group -- for future use }

  \item{se.N}{ logical for whether to estimate SE(\eqn{\hat{N}}{N-hat})
  and confidence interval }

  \item{alpha}{ alpha level for confidence intervals }

  \item{loginterval}{ logical for whether to base interval on log(N)}

  \item{keep.region}{ logical for whether to save the raster region}

  \item{nlowerbound}{ logical for whether to use n as lower bound when
    computing log interval for realised N}

  \item{RN.method}{ character string for method used to calculate
    realised N (RN) and its sampling variance. `poisson' or `MSPE'.}

  \item{pooled.RN}{ logical; if TRUE the estimate of realised N for a multi-session
  model is computed as if for combined sampling with all detectors (see Details)}
  
  \item{ncores}{ integer number of threads to be used for parallel processing}

  \item{...}{ other arguments (not used)}

}

\details{
  
  If the density surface of the fitted model is flat
  (i.e. \code{object$model$D == ~1} or \code{object$CL == TRUE}) then
  \eqn{E(N)} is simply the density multiplied by the area of \code{region},
  and the standard error is also a simple product. In the conditional
  likelihood case, the density and standard error are obtained by first
  calling \code{derived}.

  If, on the other hand, the density has been modelled then the density
  surface is predicted at each point in \code{region} and \eqn{E(N)} is
  obtained by discrete summation. Pixel size may have a minor effect on
  the result - check by varying \code{spacing}. Sampling variance is
  determined by the delta method, using a numerical approximation to the
  gradient of \eqn{E(N)} with respect to each beta parameter.

  The region may be defined as a mask object (if omitted, the mask
  component of \code{object} will be used). Alternatively, \code{region}
  may be a SpatialPolygonsDataFrame object (see package \pkg{sp}), and a
  raster mask will be constructed on the fly using the specified
  spacing. See \code{\link{make.mask}} for an example importing a
  shapefile to a SpatialPolygonsDataFrame.

  Note: The option of specifying a polygon rather than a mask for
  \code{region} does not work if the density model in \code{object} uses
  spatial covariates: these must be passed in a mask.
  
  Group-specific N has yet to be implemented.

  Population size is adjusted automatically for the number of clusters
  in `mashed' models (see \code{\link{mash}}). However, the population
  size reported is that associated with a single cluster unless
  \code{regionmask} is specified.

\code{pooled.RN = TRUE} handles the special case of a multi-session
model in which the region of interest spans several patches (i.e.,
sampling in each session is localised within \code{region}. This is not
yet fully implemented.
  
Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

Use \code{\link{par.region.N}} to apply \code{region.N} in parallel to
several models.

}

\value{

  If \code{se.N = FALSE}, the numeric value of expected population size,
  otherwise, a dataframe with rows `E.N' and `R.N', and columns as
  below.
  
\tabular{ll}{
estimate    \tab estimate of N (expected or realised, depending on row)\cr
SE.estimate \tab standard error of estimated N \cr
lcl         \tab lower 100(1--alpha)\% confidence limit \cr
ucl         \tab upper 100(1--alpha)\% confidence limit \cr
n           \tab total number of individuals detected \cr
}

For multiple sessions, the value is a list with one component per
session, each component as above.

If \code{keep.region = TRUE} then the mask object for the region is
saved as the attribute `region' (see Examples).

The area in hectares of the region is saved as attribute `regionarea'.

}

\note{

  The estimates of expected and realised \eqn{N} are generally very
  similar, or identical, but realised \eqn{N} usually has lower
  estimated variance, especially if the \eqn{n} detected animals
  comprise a large fraction. 
 
  Realised \eqn{N} is given by \eqn{R(N) = n + \int_B (1 - p.(X))D(X)
  dX}{R(N) = n + integral_over_B (1 - p.(X))D(X) dX} (the second term
  represents undetected animals). This definition strictly holds only
  when region B is at least as large as the region of integration used
  to fit the model; only with this condition can we be sure all \eqn{n}
  detected animals have centres within B. The sampling variance of
  \eqn{R(N)}, technically a mean square prediction error (Johnson et al.
  2010), is approximated by summing the expected Poisson variance of the
  true number of undetected animals and a delta-method estimate of its
  sampling variance, obtained as for \eqn{E(N)}.

  By default, a shortcut is used to compute the sampling variance of
  realised \eqn{N}. With this option (RN.method = `poisson') the
  sampling variance is the sampling variance of \eqn{E(N)} minus the
  estimate of \eqn{E(N)} (representing Poisson process variance). This
  has been found to give reliable confidence intervals in simulations
  (Efford and Fewster 2013).

  If RN.method is neither `MSPE' nor `poisson' (ignoring case) then
  the estimate of expected \eqn{N} is also used for realised \eqn{N},
  and the `poisson' shortcut variance is used.
  
  Johnson et al. (2010) use the notation \eqn{\mu(B)} for expected
  \eqn{N} and \eqn{N(B)} for realised \eqn{N} in region \eqn{B}.

  In our case, the relative SE (CV) of \eqn{\mu(B)} is the same as that
  for the estimated density \eqn{D} if \eqn{D} has been estimated using
  the Poisson distribution option in \code{secr.fit} or
  \code{derived()}. If \eqn{D} has been estimated with the binomial
  distribution option, its relative SE for simple models will be the
  same as that of \eqn{N(B)}, assuming that \eqn{B} is the full extent
  of the original mask.

}
  
\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. and Fewster, R. M. (2013) Estimating population size by
spatially explicit capture--recapture. \emph{Oikos} \bold{122}, 918--928.

Johnson, D. S., Laake, J. L. and Ver Hoef, J. M. (2010) A model-based
approach for making ecological inference from distance sampling
data. \emph{Biometrics} \bold{66}, 310--318.

}

\seealso{
  
  \code{\link{secr.fit}}, \code{\link{derived}}, \code{\link{make.mask}},
  \code{\link{expected.n}}, \code{\link{closedN}}
  
}

\examples{

\dontrun{

## routine examples using arbitrary mask from model fit
region.N(secrdemo.0)
region.N(secrdemo.CL)
region.N(ovenbird.model.D)

## region defined as vector polygon
## retain and plot region mask
temp <- region.N(possum.model.0, possumarea, spacing = 40,
    keep.region = TRUE)
temp
plot (attr(temp, "region"), type = "l")
}

}

\keyword{ models }
%E.n         \tab expected number of individuals detected from region \cr
