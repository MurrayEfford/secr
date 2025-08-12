\name{derivedDsurface}
\alias{derivedDsurface}
\alias{derivedDcoef}
\alias{derivedDcoef.secr}
\alias{derivedDcoef.secrlist}

\title{Derived Absolute Density}

\description{

A relative density model (one fitted by maximising the conditional likelihood) has 
a zero intercept for density on the default log scale (the intercept is 1.0 for the identity link). The constant \emph{k} that relates absolute density to relative density \eqn{D^\prime = k^{-1} D} is obtained as described in Efford (2025).

These functions infer the value of \emph{k} and use this to construct the predicted density surface (\code{\link{predictDsurface}} is the equivalent function for full-likelihood models).

The \code{\link{derived}} method for `secr' objects performs the same calculation as \code{derivedDcoef} and also back-transforms the intercept and computes a delta-method estimate of its variance and confidence limits. However that variance estimate is unreliable.

}

\usage{

\method{derivedDcoef}{secr}(object, se = FALSE, ...)
\method{derivedDcoef}{secrlist}(object, se = FALSE, ...)

derivedDsurface(object, mask = NULL, sessnum = NULL, groups = NULL)

}

\arguments{

  \item{object}{fitted secr relative density model}
  \item{se}{logical; if TRUE the variance of beta0 is estimated by the delta method}
  \item{mask}{new mask for which to compute densities}
  \item{sessnum}{integer session number}
  \item{groups}{character vector of covariate names to define groups (optional)}
  \item{\dots}{other arguments (not used)}

}

\details{

The theory is provided by Efford in prep. 

\code{derivedDcoef(object, se = TRUE)} is equivalent to \code{derived(object, Dweight = TRUE)}, 
although \code{derivedDcoef} returns estimates on the link scale.

}

\value{

For derivedDcoef -- 

A dataframe mimicking \code{coef(object)} with an initial row for the derived density intercept  \eqn{\beta_0}.
For an identity link, subsequent density coefficients are scaled by the derived \eqn{\beta_0}.

For derivedDsurface --

Dsurface object with class c("Dsurface", "mask", "data.frame"). 
Multi-session if object is multi-session and sessnum = NULL.

If groups are defined the result is a list of Dsurfaces.

}

\references{

Efford, M. G. In prep. SECR models for relative density. 

}

\seealso{

\code{\link{derived}}

\code{\link{predictDsurface}}

}