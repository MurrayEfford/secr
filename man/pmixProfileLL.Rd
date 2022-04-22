\name{pmixProfileLL}
\alias{pmixProfileLL}

\title{Mixture Model Check}

\description{
  
Compute the profile likelihood of a finite mixture model for a
user-specified range of values for the mixing parameter. This provides a
check on multimodality.

}

\usage{

pmixProfileLL(CH, model = list(g0 ~ h2, sigma ~ h2), CL = TRUE, pmvals = seq(0.01,
 0.99, 0.01), pmi = 5, ...)

}
\arguments{
  \item{CH}{capthist object}
  \item{model}{model as in \code{\link{secr.fit}}}
  \item{CL}{logical as in in \code{\link{secr.fit}}}
  \item{pmvals}{numeric vector of values for mixing parameter `pmix'}
  \item{pmi}{integer index of `pmix' in vector of coefficients (beta
    parameters) for the specified model}
  \item{\dots}{other arguments passed to \code{\link{secr.fit}}}
}

\details{
  
  See \href{https://www.otago.ac.nz/density/pdfs/secr-finitemixtures.pdf}{secr-finitemixtures.pdf}.

  Choosing the wrong value for pmi results in the error message "invalid fixed beta - require NP-vector". The easiest way to find the value of \code{pmi} is to inspect the
  output from a previously fitted mixture model - either count the coefficients
  or check fit$parindx$pmix (for a model named `fit'). It is assumed that `pmix' is the last real
  parameter in the model, and that pmix is constant.
  
}

\value{

  Numeric vector of profile likelihoods.
  
}

\note{

  This is slow to execute and the results are hard to interpret. Use
  only if you are confident.

}

\examples{

\dontrun{

pmvals <- seq(0.02,0.99,0.02)
mask <- make.mask(traps(ovenCH[[1]]), nx = 32, buffer = 100)

## only g0 ~ h2, so reduce pmi from 5 to 4
outPL <- pmixProfileLL(ovenCH[[1]], model = list(g0~h2), 
    mask = mask, pmvals, CL = TRUE, trace = FALSE, pmi = 4) 
    
plot(pmvals, outPL, xlim = c(0,1),
xlab = 'Fixed pmix', ylab = 'Profile log-likelihood')

}
}

\keyword{ model }
