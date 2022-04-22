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

# Deprecated in 4.5.4

ip.secr (capthist, predictorfn = pfn, predictortype = "null", detectfn = 0,
    mask = NULL, start = NULL, boxsize = 0.2, boxsize2 = boxsize, centre = 3,  
    min.nsim = 10, max.nsim = 2000, CVmax = 0.002, var.nsim = 1000, maxbox = 5, 
    maxtries = 2, ncores = 1, seed = NULL, trace = TRUE, ...)

pfn(capthist, N.estimator = c("n", "null","zippin","jackknife") )

\code{model.average} is replaced by a method for `secr' and `secrlist' objects of the generic \code{\link{modelAverage}}. The internal code is essentially the same for \code{model.average} and \code{modelAverage.secrlist}. The generic avoids a name conflict with RMark and is also used in \pkg{openCR}.
 
\code{ip.secr} will be replaced by function \code{ipsecr.fit} in the forthcoming package \pkg{ipsecr}. 
 
}

\seealso{

\code{\link{secr-defunct}}, 
\code{\link{modelAverage}}, 
\code{\link{model.average-deprecated}}, 
\code{\link{ip.secr-deprecated}}


}

