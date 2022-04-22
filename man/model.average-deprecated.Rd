\name{model.average-deprecated}
\alias{model.average-deprecated}
\alias{model.average}

\title{Old Model Average Function}

\description{

Deprecated in 4.5.0 and will be removed in future versions of \pkg{secr}.

}

\usage{

model.average (..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, dmax = 10, covar = FALSE, average = c('link', 'real'), 
    criterion = c('AICc','AIC'), CImethod = c('Wald', 'MATA'))

}

\arguments{
  \item{\dots}{ \code{secr} or \code{secrlist} objects }
  \item{realnames}{ character vector of real parameter names }
  \item{betanames}{ character vector of beta parameter names }
  \item{newdata}{ optional dataframe of values at which to evaluate models }
  \item{alpha}{ alpha level for confidence intervals }
  \item{dmax}{ numeric, the maximum AIC or AICc difference for inclusion in confidence set }
  \item{covar}{ logical, if TRUE then return variance-covariance matrix }
  \item{average}{ character string for scale on which to average real
    parameters }
  \item{criterion}{ character, information criterion to use for model weights}  
  \item{CImethod}{ character, type of confidence interval (see Details)}
}

\details{

\code{model.average} is replaced by a method for `secr' and `secrlist' objects of the generic \code{\link{modelAverage}}. The internal code is essentially the same for \code{model.average} and \code{modelAverage.secrlist}. The generic avoids a name conflict with RMark and is also used in \pkg{openCR}.
 
}

\seealso{

\code{\link{secr-deprecated}}

}

