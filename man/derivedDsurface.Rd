\name{derivedDsurface}
\alias{derivedDsurface}
\alias{derivedDbeta0}

\title{Derived Absolute Density}

\description{

A relative density model (one fitted by maximising the conditional likelihood) has 
a zero intercept for density on the default log scale (the intercept is 1.0 for the identity link). The constant \emph{k} that relates absolute density to relative density \eqn{D^\prime = k^{-1} D} is obtained as described in Efford (2025).

These functions infer the value of \emph{k} and use this to construct the predicted density surface (\code{\link{predictDsurface}} is the equivalent function for full-likelihood models).

The \code{\link{derived}} method for `secr' objects performs the same calculation as \code{derivedDbeta0} and also back-transforms the intercept and computes a delta-method estimate of its variance and confidence limits. \code{derivedDbeta0} is mostly for internal use.

}

\usage{

derivedDbeta0(object, sessnum = 1, groups = NULL, se.beta0 = FALSE)

derivedDsurface(object, mask = NULL, sessnum = NULL, groups = NULL)

}

\arguments{

  \item{object}{fitted secr relative density model}
  \item{mask}{new mask for which to compute densities}
  \item{sessnum}{integer session number}
  \item{groups}{character vector of covariate names to define groups (optional)}
  \item{se.beta0}{logical; if TRUE the variance is estimated by the delta method}

}

\details{

The theory is provided by Efford in prep. 

\code{derivedDbeta0(object, se.beta0 = TRUE)} is equivalent to \code{derived(object, Dweight = TRUE)}, although \code{derivedDbeta0} returns estimates on the link scale.

}

\value{

For derivedDbeta0 -- 

The derived intercept (and optionally the SE) of the density model on the link scale i.e. \eqn{f(k)} where \eqn{f} is the link function.

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