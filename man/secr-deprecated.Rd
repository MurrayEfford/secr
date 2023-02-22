\name{secr-deprecated}
\alias{secr-deprecated}

\title{Deprecated Functions in Package \pkg{secr}}

\description{

These functions will be removed from future versions of \pkg{secr}.

}

\details{

# Deprecated in 4.5.0

model.average (..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, dmax = 10, covar = FALSE, average = c('link', 'real'), 
    criterion = c('AICc','AIC'), CImethod = c('Wald', 'MATA'))

\code{model.average} is replaced by a method for `secr' and `secrlist' objects of the generic \code{\link{modelAverage}}. The internal code is essentially the same for \code{model.average} and \code{modelAverage.secrlist}. The generic avoids a name conflict with RMark and is also used in \pkg{openCR}.
 
}

\seealso{

\code{\link{secr-defunct}}, 
\code{\link{modelAverage}}, 
\code{\link{model.average-deprecated}}, 


}

