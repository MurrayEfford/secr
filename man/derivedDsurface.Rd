\name{derivedDsurface}
\alias{derivedDsurface}
\alias{derivedIntercept}

\title{Derived Absolute Density}

\description{

A relative density model (one fitted by maximising the conditional likelihood) has 
a zero intercept for density on the default log scale (the intercept is 1.0 for the identity link). The constant \emph{k} that relates absolute density to relative density \eqn{D^\prime = k^{-1} D} is obtained as described in Efford (2025).

These functions infer the value of \emph{k} and use this in a call to \code{\link{predictDsurface}}.

}

\usage{

derivedIntercept(object, sessnum = 1, groups = NULL, Dweight = TRUE)

derivedDsurface(object, mask = NULL, sessnum = NULL, groups = NULL)

}

\arguments{

  \item{object}{fitted secr relative density model}
  \item{mask}{new mask for which to compute densities}
  \item{sessnum}{integer session number}
  \item{groups}{character vector of covariate names to define groups (optional)}
  \item{Dweight}{logical; FALSE for testing}

}

\details{

The theory is provided by Efford in prep. 

\code{derivedIntercept} is defined for detection models lacking individual covariates.

}

\value{

For derivedIntercept -- 

The derived intercept of the density model on the link scale i.e. \eqn{f(k)} where \eqn{f} is the link function.

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