\name{CV}
\alias{CV}
\alias{CVa}
\alias{CVa0}

\title{ Coefficient of Variation }

\description{

  The coefficient of variation of effective sampling area predicts the
  bias in estimated density (Efford and Mowat 2014). These functions
  assist its calculation from fitted finite mixture models.
  
}

\usage{

CV(x, p, na.rm = FALSE)
CVa0(object, ...)
CVa(object, sessnum = 1, ...)

}

\arguments{
  
  \item{x}{vector of numeric values }
  \item{p}{vector of class probabilities}
  \item{na.rm}{logical; if TRUE missing values are dropped from x}
  \item{object}{fitted secr finite mixture model}
  \item{sessnum}{integer sequence number of session to analyse }
  \item{\dots}{other arguments passed to predict.secr (e.g.,
    \code{newdata})}
}
  
\details{

    \code{CV} computes the coefficient of variation of \code{x}. If
    \code{p} is provided then the distribution is assumed to be
    discrete, with support \code{x} and class membership probabilities
    \code{p} (scaled automatically to sum to 1.0).

    \code{CVa} computes CV(\eqn{a}) where \eqn{a} is the effective
    sampling area of Borchers and Efford (2008).
    
    \code{CVa0} computes CV(a0) where a0 is the single-detector sampling
    area defined as \eqn{a_0 = 2 \pi \lambda_0 \sigma^2}{a0 = 2 \pi
    lambda0.sigma^2} (Efford and Mowat 2014); a0 is a convenient
    surrogate for \emph{a}, the effective sampling area. CV(a0) uses
    either the fitted MLE of a0 (if the a0 parameterization has been
    used), or a0 computed from the estimates of lambda0 and sigma.

    \code{CVa} and \code{CVa0} do not work for models with individual
    covariates.

}
\value{

  Numeric

}

\note{

  Do not confuse the function CVa with the estimated relative standard
  error of the estimate of a from \code{\link{derived}}, also labelled CVa
  in the output. The relative standard error RSE is often labelled CV
  in the literature on capture--recapture, but this can cause unnecessary 
  confusion. See also \code{\link{RSE}}.

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. and Mowat, G. (2014) Compensatory heterogeneity in
capture--recapture data. \emph{Ecology} \bold{95}, 1341--1348.  

}

\seealso{

  \code{\link{CVpdot}}, 
  \code{\link{derived}}, 
  \code{\link{details}}, 
  \code{\link{RSE}}
  
}

\examples{

\dontrun{

## housemouse model
morning <- subset(housemouse, occ = c(1,3,5,7,9))
msk <- make.mask((traps(morning)), nx = 32) 
morning.h2   <- secr.fit(morning, buffer = 20, model = list(g0~h2), mask = msk, 
    trace = FALSE)
CVa0(morning.h2 )

}

}

\keyword{ manip }