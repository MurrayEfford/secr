\name{derivedDsurface}
\alias{derivedDsurface}

\title{Derived Absolute Density}

\description{

A relative density model (one fitted with details argument relativeD = TRUE) has 
a zero intercept for density on the log scale. This function infers the value of 
the intercept that yields absolute densities (internal function \code{derivedIntercept}) 
and uses this in a call to \code{\link{predictDsurface}}.

}

\usage{

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

Dsurface object with class c("Dsurface", "mask", "data.frame"). The computed value of the log intercept is stored as the attribute \code{derivedIntercept}.

}

\references{

Efford, M. G. In prep. SECR models for relative density. 

}

\seealso{

\code{\link{predictDsurface}}, \link{details}

}