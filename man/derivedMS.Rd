\name{derived}
\alias{derived}
\alias{esa}
\alias{esa.secr}
\alias{derived.secr}
\alias{derived.secrlist}

\title{ Derived Parameters of Fitted SECR Model }

\description{

Compute derived parameters of spatially explicit capture-recapture model. 
Density is a derived parameter when a model is fitted by maximizing the conditional likelihood. So also is the effective sampling area (in the sense of Borchers and Efford 2008).

}

\usage{

derived(object, ...)

\method{derived}{secr}(object, sessnum = NULL, groups = NULL, alpha = 0.05, 
    se.esa = FALSE, se.D = TRUE, loginterval = TRUE, distribution = NULL, 
    ncores = NULL, bycluster = FALSE, ...)

\method{derived}{secrlist}(object, sessnum = NULL, groups = NULL, alpha = 0.05, 
    se.esa = FALSE, se.D = TRUE, loginterval = TRUE, distribution = NULL, 
    ncores = NULL, bycluster = FALSE, ...)

\method{esa}{secr}(object, sessnum = 1, beta = NULL, real = NULL, noccasions = NULL, 
    ncores = NULL, ...)

}

\arguments{
  \item{object}{ \code{secr} object output from \code{secr.fit}, or an
  object of class \code{c("secrlist", "list")} }
  \item{sessnum}{ index of session in object$capthist for which output required }
  \item{groups}{ vector of covariate names to define group(s) (see Details) }
  \item{alpha}{ alpha level for confidence intervals }
  \item{se.esa}{ logical for whether to calculate SE(mean(esa)) }
  \item{se.D}{ logical for whether to calculate SE(D-hat) }
  \item{loginterval}{ logical for whether to base interval on log(D)}
  \item{distribution}{ character string for distribution of the number of individuals detected }
  \item{ncores}{integer number of threads used for parallel processing}
  \item{bycluster}{logical; if TRUE results are reported separately for each cluster of detectors}
  \item{beta}{ vector of fitted parameters on transformed (link) scale }
  \item{real}{ vector of `real' parameters } 
  \item{noccasions}{ integer number of sampling occasions (see Details) } 
  \item{...}{other arguments passed to \code{\link[nlme]{fdHess}} or ignored (esa)}
}

\details{

The derived estimate of density is a Horvitz-Thompson-like estimate:
\deqn{\hat{D} = \sum\limits _{i=1}^{n} {a_i (\hat{\theta})^{-1}}}{ D-hat
= sum( a_i (theta-hat)^--1)} where \eqn{a_i (\hat{\theta})}{% a_i
(theta-hat)} is the estimate of effective sampling area for animal
\eqn{i} with detection parameter vector \eqn{\theta}.

A non-null value of the argument \code{distribution} overrides the value
in \code{object$details}. The sampling variance of \eqn{\hat{D}}{D-hat}
from \code{secr.fit} by default is spatially unconditional
(\code{distribution = "Poisson"}). For sampling variance conditional on the population of the
habitat mask (and therefore dependent on the mask area), specify
\code{distribution = "binomial"}. The equation for the conditional
variance includes a factor \eqn{(1 - a/A)} that disappears in the
unconditional (Poisson) variance (Borchers and Efford 2007). Thus the
conditional variance is always less than the unconditional variance. The
unconditional variance may in turn be an overestimate or (more likely)
an underestimate if the true spatial variance is non-Poisson.

Derived parameters may be estimated for population subclasses (groups)
defined by the user with the \code{groups} argument. Each named factor
in \code{groups} should appear in the covariates dataframe of
object$capthist (or each of its components, in the case of a
multi-session dataset).

\code{esa} is used by \code{derived} to compute individual-specific
effective sampling areas: \deqn{a_i (\hat{\theta}) = \int _A \:
p.(\mathbf{X};\mathbf{z}_i, \mathbf{\hat{\theta}}) \; \mathrm{d}
\mathbf{X}}{a_i = integral p.(X; z_i, theta_i) dX} where
\eqn{p.(\mathbf{X})}{p.(X)} is the probability an individual at X is
detected at least once and the \eqn{\mathbf{z}_i}{z_i} are optional
individual covariates. Integration is over the area \eqn{A} of the
habitat mask.

The argument \code{noccasions} may be used to vary the number of
sampling occasions; it works only when detection parameters are constant
across individuals and across time.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

The effective sampling area `esa' (\eqn{a(\hat{\theta})}{a(theta-hat)})
reported by \code{derived} is equal to the harmonic mean of the
\eqn{a_i (\hat{\theta})}{a_i (theta-hat)} (arithmetic
mean prior to version 1.5). The sampling variance of
\eqn{a(\hat{\theta})}{a(theta-hat)} is estimated by
\deqn{\widehat{\mbox{var}}(a(\hat{\theta})) = \hat{G}_\theta^T
\hat{V}_\theta \hat{G}_\theta, }{ var(a(theta)) = G-hat^T V-hat G-hat,}
where \eqn{\hat{V}}{V-hat} is the asymptotic estimate of the
variance-covariance matrix of the beta detection parameters
(\eqn{\theta}{theta}) and \eqn{\hat{G}}{G-hat} is a numerical estimate
of the gradient of \eqn{a(\theta)}{a(theta)} with respect to
\eqn{\theta}{theta}, evaluated at \eqn{\hat{\theta}}{theta-hat}.

A 100(1--alpha)\% asymptotic confidence interval is reported for
density. By default, this is asymmetric about the estimate because the
variance is computed by backtransforming from the log scale. You may 
also choose a symmetric interval (variance obtained on natural scale). 

The vector of detection parameters for \code{esa} may be specified via
\code{beta} or \code{real}, with the former taking precedence. If
neither is provided then the fitted values in \code{object$fit$par} are
used. Specifying \code{real} parameter values bypasses the various
linear predictors. Strictly, the `real' parameters are for a naive
capture (animal not detected previously).

The computation of sampling variances is relatively slow and may be
suppressed with \code{se.esa} and \code{se.D} as desired.

For computing \code{derived} across multiple models in parallel see
\code{\link{par.derived}}.

From \pkg{secr} 5.0.0 the \dots argument may be used to control the step size 
(.relStep) used by \code{\link[nlme]{fdHess}} when estimating gradients for SE(D) and SE(esa).

}

\value{

Dataframe with one row for each derived parameter (`esa', `D') and
columns as below

\tabular{ll}{
estimate    \tab estimate of derived parameter \cr
SE.estimate \tab standard error of the estimate \cr
lcl         \tab lower 100(1--alpha)\% confidence limit \cr
ucl         \tab upper 100(1--alpha)\% confidence limit \cr
CVn         \tab relative SE of number observed (Poisson or binomial assumption) \cr
CVa         \tab relative SE of effective sampling area \cr
CVD         \tab relative SE of density estimate \cr
}

For a multi-session or multi-group analysis the value is a list
with one component for each session and group.

The result will also be a list if \code{object} is an `secrlist'.

}

\references{
Borchers, D. L. and Efford, M. G. (2007) Supplements to Biometrics paper. Available online at \url{https://www.otago.ac.nz/density/}.

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture--recapture studies. \emph{Biometrics}, \bold{64}, 377--385.

Huggins, R. M. (1989) On the statistical analysis of capture experiments. \emph{Biometrika} \bold{76}, 133--140.

}

\seealso{

  \code{\link{predict.secr}},
  \code{\link{print.secr}},
  \code{\link{secr.fit}},
  \code{\link{empirical.varD}}

}

\note{

Before version 2.1, the output table had columns for `varcomp1' (the variance in \eqn{\hat{D}}{D-hat} due to variation in \eqn{n}, i.e.,
Huggins' \eqn{s^2}), and `varcomp2' (the variance in \eqn{\hat{D}}{D-hat} due to uncertainty in estimates of detection parameters).

These quantities are related to CVn and CVa as follows:

\eqn{\mbox{CVn} = \sqrt{ \mbox{varcomp1} } / \hat{D} }{ CVn = sqrt(varcomp1) / D-hat}

\eqn{\mbox{CVa} = \sqrt{ \mbox{varcomp2} } / \hat{D} }{ CVa = sqrt(varcomp2) / D-hat}

}

\section{Warning}{
%2014-06-08
\code{derived()} may be applied to detection models fitted by maximizing the full likelihood (\code{CL = FALSE}). However, models using g (groups) will not be handled correctly.

}

\examples{

\dontrun{
## extract derived parameters from a model fitted previously
## by maximizing the conditional likelihood 
derived (secrdemo.CL)

## what happens when sampling variance is conditional on mask N?
derived(secrdemo.CL, distribution = "binomial")
## fitted g0, sigma
esa(secrdemo.CL)
## force different g0, sigma
esa(secrdemo.CL, real = c(0.2, 25))
}

}

\keyword{ models }

%varcomp1    \tab variance due to variation in \eqn{n} (Huggins' \eqn{s^2}) \cr
%varcomp2    \tab variance due to uncertainty in estimates of detection parameters \cr
