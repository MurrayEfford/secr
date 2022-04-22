\name{empirical.varD}
\alias{empirical.varD}
\alias{derivednj}
\alias{derivedSession}
\alias{derivedMash}
\alias{derivedCluster}
\alias{derivedExternal}
\alias{derivedSystematic}

\title{

Empirical Variance of H-T Density Estimate

}

\description{

  Compute Horvitz-Thompson-like estimate of population density from a
  previously fitted spatial detection model, and estimate its sampling
  variance using the empirical spatial variance of the number observed
  in replicate sampling units. Wrapper functions are provided for
  several different scenarios, but all ultimately call
  \code{derivednj}. The function \code{derived} also computes
  Horvitz-Thompson-like estimates, but it assumes a Poisson or binomial
  distribution of total number when computing the sampling variance.

}

\usage{

derivednj ( nj, esa, se.esa = NULL, method = c("SRS", "R2", "R3", "local",
    "poisson", "binomial"), xy = NULL, alpha = 0.05, loginterval = TRUE, 
    area = NULL, independent.esa = FALSE )

derivedMash ( object, sessnum = NULL,  method = c("SRS", "local"),
    alpha = 0.05, loginterval = TRUE)

derivedCluster ( object, method = c("SRS", "R2", "R3", "local", "poisson", "binomial"),
    alpha = 0.05, loginterval = TRUE)

derivedSession ( object,  method = c("SRS", "R2", "R3", "local", "poisson", "binomial"), 
    xy = NULL, alpha = 0.05, loginterval = TRUE, area = NULL, independent.esa = FALSE )

derivedExternal ( object, sessnum = NULL, nj, cluster, buffer = 100,
    mask = NULL, noccasions = NULL,  method = c("SRS", "local"), xy = NULL,
    alpha = 0.05, loginterval = TRUE)

derivedSystematic( object, xy, design = list(), basenx = 10, df = 9, extrapolate = TRUE,
    alpha = 0.05, loginterval = TRUE, independent.esa = FALSE, keep = FALSE,
    ncores = NULL)

}

\arguments{
  \item{object}{fitted secr model }
  \item{nj}{ vector of number observed in each sampling unit (cluster)}
  \item{esa}{estimate of effective sampling area (\eqn{\hat{a}}{a-hat})}
  \item{se.esa}{ estimated standard error of effective sampling
    area (\eqn{\widehat{SE}(\hat{a})}{SE-hat(a-hat)})}
  \item{method}{ character string `SRS' or `local' }
  \item{xy}{ dataframe of x- and y- coordinates (\code{method = "local"} only)}
  \item{alpha}{ alpha level for confidence intervals }
  \item{loginterval}{ logical for whether to base interval on log(N)}
  \item{area}{ area of region for method = "binomial" (hectares) }
  \item{independent.esa}{logical; controls variance contribution from esa (see Details) }
  \item{sessnum}{ index of session in object$capthist for which output required }
  \item{cluster}{ `traps' object for a single cluster}
  \item{buffer}{ width of buffer in metres (ignored if \code{mask}
  provided) } 
  \item{mask}{ mask object for a single cluster of detectors }
  \item{noccasions}{number of occasions (for \code{nj}) }
  
  \item{design}{list specifying systematic design (see Details)}
  \item{basenx}{integer number of basis grid points in x-dimension}
  \item{df}{integer number of degrees of freedom for gam}
  \item{extrapolate}{ logical; if FALSE then boxlet p values are inferred from nearest point inside convex hull of grid}
  \item{keep}{logical; if TRUE then derivedSystematic saves key intermediate values as attributes}
  \item{ncores}{integer}
}

\details{

  \code{derivednj} accepts a vector of counts (\code{nj}), along with
  \eqn{\hat{a}}{a-hat} and \eqn{\widehat{SE}(\hat{a})}{SE(a-hat)}.  The
  argument \code{esa} may be a scalar or (if se.esa is NULL) 
  a 2-column matrix with \eqn{\hat{a_j}}{a_j-hat} and
  \eqn{\widehat{SE}(\hat{a_j})}{SE(a_j-hat)} for each replicate \eqn{j} (row). 
  In the special case that \code{nj}
  is of length 1, or \code{method} takes the values `poisson' or
  `binomial', the variance is computed using a theoretical variance
  rather than an empirical estimate. The value of \code{method}
  corresponds to `distribution' in \code{\link{derived}}, and defaults to
  `poisson'. For \code{method = 'binomial'} you must specify \code{area}
  (see Examples).

  If \code{independent.esa} is TRUE then independence is assumed among 
  cluster-specific estimates of esa, and esa variances are summed. The default 
  is a weighted sum leading to higher overall variance.

  \code{derivedCluster} accepts a model fitted to data from clustered
  detectors; each \link{cluster} is interpreted as a replicate
  sample. It is assumed that the sets of individuals sampled by
  different clusters do not intersect, and that all clusters have the
  same geometry (spacing, detector number etc.).
  
  \code{derivedMash} accepts a model fitted to clustered data that have
  been `mashed' for fast processing (see \code{\link{mash}}); each
  cluster is a replicate sample: the function uses the vector of cluster
  frequencies (\eqn{n_j}) stored as an attribute of the mashed
  \code{capthist} by \code{mash}.

  \code{derivedExternal} combines detection parameter estimates from a
  fitted model with a vector of frequencies \code{nj} from replicate
  sampling units configured as in \code{cluster}. Detectors in
  \code{cluster} are assumed to match those in the fitted model with
  respect to type and efficiency, but sampling duration
  (\code{noccasions}), spacing etc. may differ. The \code{mask} should
  match \code{cluster}; if \code{mask} is missing, one will be
  constructed using the \code{buffer} argument and defaults from
  \code{make.mask}.

  \code{derivedSession} accepts a single fitted model that must span
  multiple sessions; each session is interpreted as a replicate sample.

  Spatial variance is calculated by one of these methods
  
  \tabular{llll}{
  Method \tab Description \cr

\code{"SRS"} \tab simple random sampling with identical clusters \cr
\code{"R2"} \tab variable cluster size cf Thompson (2002:70) estimator for line transects \cr
\code{"R3"} \tab variable cluster size cf Buckland et al. (2001)\cr
\code{"local"} \tab neighbourhood variance estimator (Stevens and Olsen 2003) SUSPENDED in 4.4.7 \cr 
\code{"poisson"} \tab theoretical (model-based) variance \cr
\code{"binomial"} \tab theoretical (model-based) variance in given \code{area}\cr
}

The weighted options R2 and R3 substitute \eqn{\hat{a_j}}{a_j-hat} for line length \eqn{l_k}{l_k} in the corresponding formulae of Fewster et al. (2009, Eq 3,5). Density is estimated by \eqn{D = n/A}{\hat D = n/A} where \eqn{A = \sum a_j}{A = \sum a_j}. The variance of \eqn{A} is estimated as the sum of the cluster-specific variances, assuming independence among clusters. Fewster et al. (2009) found that an alternative estimator for line transects derived by Thompson (2002) performed better when there were strong density gradients correlated with line length (R2 in Fewster et al. 2009, Eq 3).

The neighborhood variance estimator is implemented in package \pkg{spsurvey} and was originally proposed for generalized random tessellation stratified (GRTS) samples.  For `local' variance
  estimates, the centre of each replicate must be provided in \code{xy},
  except where centres may be inferred from the data. It is unclear whether `local' can or should be used when clusters vary in size.

\code{derivedSystematic} implements the 'boxlet' variance estimator of Fewster (2011) for systematic designs using clustered detectors (an alternative to \code{derivedCluster} and \code{derivedSessions}). The method is experimental in secr 3.2.0 and may change. The `design' argument is a list with components corresponding to arguments of \code{\link{make.systematic}}, (\code{n} and \code{origin} are ignored if provided):

\tabular{llll}{
Component \tab Description \cr
\code{cluster} \tab traps object for a single cluster \cr
\code{region} \tab 2-column matrix or SpatialPolygons \cr
\code{spacing} \tab spacing between cluster origins \cr
\code{...} \tab other arguments passed to \code{\link{trap.builder}} \cr
\tab e.g. \code{edgemethod}, \code{exclude}, \code{exclmethod} \cr
}

If \code{region} is omitted from \code{design} then an attempt will be made to retrieve it from the mask attribute of \code{object} (this works if the call to \code{\link{make.mask}} used \code{keep.poly = TRUE}).

}

\value{

Dataframe with one row for each derived parameter (`esa', `D') and
columns as below

\tabular{ll}{
estimate    \tab estimate of derived parameter \cr
SE.estimate \tab standard error of the estimate \cr
lcl         \tab lower 100(1--alpha)\% confidence limit \cr
ucl         \tab upper 100(1--alpha)\% confidence limit \cr
CVn         \tab relative SE of number observed (across sampling units) \cr
CVa         \tab relative SE of effective sampling area \cr
CVD         \tab relative SE of density estimate \cr
}

}

\references{

Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L., Borchers,
D. L. and Thomas, L. (2001) \emph{Introduction to Distance Sampling:
Estimating Abundance of Biological Populations}. Oxford University
Press, Oxford.

Fewster, R. M. (2011) Variance estimation for systematic designs in spatial 
  surveys. \emph{Biometrics} \bold{67}, 1518--1531.

Fewster, R. M., Buckland, S. T., Burnham, K. P., Borchers, D. L., Jupp, P. E., 
Laake, J. L. and Thomas, L. (2009) Estimating the encounter rate variance in 
distance sampling. \emph{Biometrics} \bold{65}, 225--236. 

Stevens, D. L. Jr and Olsen, A. R. (2003) Variance estimation for
spatially balanced samples of environmental resources.
\emph{Environmetrics} \bold{14}, 593--610.

Thompson, S. K. (2002) \emph{Sampling}. 2nd edition. Wiley, New York.

}

\seealso{
\code{\link{derived}}, \code{\link{esa}} 
}

\note{

  The variance of a Horvitz-Thompson-like estimate of density may be
  estimated as the sum of two components, one due to uncertainty in the
  estimate of effective sampling area (\eqn{\hat{a}}{a-hat}) and the
  other due to spatial variance in the total number of animals \eqn{n}
  observed on \eqn{J} replicate sampling units (\eqn{n =
    \sum_{j=1}^{J}{n_j}}{sum(n_j)}). We use a delta-method approximation
  that assumes independence of the components:

  \deqn{ \widehat{\mbox{var}}(\hat{D}) = \hat{D}^2
    \{\frac{\widehat{\mbox{var}}(n)}{n^2} +
    \frac{\widehat{\mbox{var}}(\hat{a})}{\hat{a}^2}\} }{var(D-hat) =
    D-hat^2 (var(n)/n^2 + var(a-hat)/a-hat^2) }

where \eqn{\widehat{\mbox{var}}(n) = \frac{J}{J-1}
  \sum_{j=1}^{J}(n_j-n/J)^2}{var(n) = J/(J-1).sum((n_j - n/J)^2)}. The
  estimate of \eqn{\mbox{var}(\hat{a})}{var(a-hat)} is model-based while
  that of \eqn{\mbox{var}(n)}{var(n)} is design-based. This formulation follows
  that of Buckland et al. (2001, p. 78) for conventional distance
  sampling. Given sufficient independent replicates, it is a robust way
  to allow for unmodelled spatial overdispersion.

  There is a complication in SECR owing to the fact that
  \eqn{\hat{a}}{a-hat} is a derived quantity (actually an integral)
  rather than a model parameter. Its sampling variance
  \eqn{\mbox{var}(\hat{a})}{var(a-hat)} is estimated indirectly in
  \pkg{secr} by combining the asymptotic estimate of the covariance
  matrix of the fitted detection parameters \eqn{\theta}{theta} with a
  numerical estimate of the gradient of \eqn{a(\theta)}{a(theta)} with
  respect to \eqn{\theta}{theta}. This calculation is performed in
  \code{\link{derived}}.
  
}

\examples{

## The `ovensong' data are pooled from 75 replicate positions of a
## 4-microphone array. The array positions are coded as the first 4
## digits of each sound identifier. The sound data are initially in the
## object `signalCH'. We first impose a 52.5 dB signal threshold as in
## Dawson & Efford (2009, J. Appl. Ecol. 46:1201--1209). The vector nj
## includes 33 positions at which no ovenbird was heard. The first and
## second columns of `temp' hold the estimated effective sampling area
## and its standard error.

\dontrun{

signalCH.525 <- subset(signalCH, cutval = 52.5)
nonzero.counts <- table(substring(rownames(signalCH.525),1,4))
nj <- c(nonzero.counts, rep(0, 75 - length(nonzero.counts)))
temp <- derived(ovensong.model.1, se.esa = TRUE)
derivednj(nj, temp["esa",1:2])

## The result is very close to that reported by Dawson & Efford
## from a 2-D Poisson model fitted by maximizing the full likelihood.

## If nj vector has length 1, a theoretical variance is used...
msk <- ovensong.model.1$mask
A <- nrow(msk) * attr(msk, "area")
derivednj (sum(nj), temp["esa",1:2], method = "poisson")
derivednj (sum(nj), temp["esa",1:2], method = "binomial", area = A)

## Set up an array of small (4 x 4) grids,
## simulate a Poisson-distributed population,
## sample from it, plot, and fit a model.
## mash() condenses clusters to a single cluster

testregion <- data.frame(x = c(0,2000,2000,0),
    y = c(0,0,2000,2000))
t4 <- make.grid(nx = 4, ny = 4, spacing = 40)
t4.16 <- make.systematic (n = 16, cluster = t4,
    region = testregion)
popn1 <- sim.popn (D = 5, core = testregion,
    buffer = 0)
capt1 <- sim.capthist(t4.16, popn = popn1)
fit1 <- secr.fit(mash(capt1), CL = TRUE, trace = FALSE)

## Visualize sampling
tempmask <- make.mask(t4.16, spacing = 10, type =
    "clusterbuffer")
plot(tempmask)
plot(t4.16, add = TRUE)
plot(capt1, add = TRUE)

## Compare model-based and empirical variances.
## Here the answers are similar because the data
## were simulated from a Poisson distribution,
## as assumed by \code{derived}

derived(fit1)
derivedMash(fit1)

## Now simulate a patchy distribution; note the
## larger (and more credible) SE from derivedMash().

popn2 <- sim.popn (D = 5, core = testregion, buffer = 0,
    model2D = "hills", details = list(hills = c(-2,3)))
capt2 <- sim.capthist(t4.16, popn = popn2)
fit2 <- secr.fit(mash(capt2), CL = TRUE, trace = FALSE)
derived(fit2)
derivedMash(fit2)

## The detection model we have fitted may be extrapolated to
## a more fine-grained systematic sample of points, with
## detectors operated on a single occasion at each...
## Total effort 400 x 1 = 400 detector-occasions, compared
## to 256 x 5 = 1280 detector-occasions for initial survey.

t1 <- make.grid(nx = 1, ny = 1)
t1.100 <- make.systematic (cluster = t1, spacing = 100,
    region = testregion)
capt2a <- sim.capthist(t1.100, popn = popn2, noccasions = 1)
## one way to get number of animals per point
nj <- attr(mash(capt2a), "n.mash")
derivedExternal (fit2, nj = nj, cluster = t1, buffer = 100,
    noccasions = 1)

## Review plots
library(MASS)
base.plot <- function() {
    eqscplot( testregion, axes = FALSE, xlab = "",
        ylab = "", type = "n")
    polygon(testregion)
}
par(mfrow = c(1,3), xpd = TRUE, xaxs = "i", yaxs = "i")
base.plot()
plot(popn2, add = TRUE, col = "blue")
mtext(side=3, line=0.5, "Population", cex=0.8, col="black")
base.plot()
plot (capt2a, add = TRUE,title = "Extensive survey")
base.plot()
plot(capt2, add = TRUE, title = "Intensive survey")
par(mfrow = c(1,1), xpd = FALSE, xaxs = "r", yaxs = "r")  ## defaults


## Weighted variance

derivedSession(ovenbird.model.1, method = "R2")

}

}

\keyword{ models }