\name{derivedDsurface}
\alias{derivedDsurface}
\alias{derivedIntercept}

\title{Derived Absolute Density}

\description{

A relative density model (one fitted by maximising the conditional likelihood) has 
a zero intercept for density on the default log scale (the intercept is 1.0 for the identity link). The constant \emph{k} that relates absolute density to relative density is obtained as described in Efford (2025).

These functions infer the value of \emph{k} and use this in a call to \code{\link{predictDsurface}}.

}

\usage{

derivedIntercept(object, sessnum = 1)

derivedDsurface(object, mask = NULL, sessnum = 1)

}

\arguments{
  \item{object}{fitted secr relative density model}
  \item{mask}{new mask for which to compute densities}
  \item{sessnum}{integer session number}
}

\details{

The theory is provided by Efford in prep. The present implementation is limited 
with respect to models with multiple sessions and groups.

}

\value{

For derivedIntercept -- 

Numeric intercept of the derived density model on the link scale.

For derivedDsurface --

Dsurface object with class c("Dsurface", "mask", "data.frame"). 

}

\references{

Efford, M. G. In prep. SECR models for relative density. 

}

\seealso{

\code{\link{predictDsurface}}

}