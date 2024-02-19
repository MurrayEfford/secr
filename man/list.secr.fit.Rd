\name{list.secr.fit}
\alias{list.secr.fit}

\title{Fit Multiple SECR Models}

\description{
  
  This function is a wrapper for \code{\link{secr.fit}} that allows multiple models to be fitted.

}

\usage{

list.secr.fit (arglist, trace = FALSE, prefix = "", seed = NULL, save.intermediate = FALSE)

}

\arguments{
  
  \item{arglist}{list of argument lists for \code{secr.fit} or a character
    vector naming such lists}
  
  \item{seed}{integer pseudorandom number seed}

  \item{trace}{logical; if TRUE intermediate output may be logged}
  
  \item{prefix}{character prefix for names of output}
  
  \item{save.intermediate}{logical; if TRUE then each fit is saved to an external file}

}

\details{

Names are created by prefixing \code{prefix} to the names of \code{argslist}. If \code{trace} is TRUE then the total
execution time and finish time are displayed.

  \code{trace} overrides any settings in \code{arglist}.
  
  \code{save.intermediate} causes each fit to be saved to a file with extension .RData.

}

\value{

An \code{\link{secrlist}} of model fits (see \code{\link{secr.fit}}). 

}

\seealso{
  \code{\link{secr.fit}},
  \code{\link{AIC.secr}},
  \code{\link{predict.secr}},
}

\note{

This function replaces the previous function \code{par.secr.fit}: since the introduction of multi-threading in \pkg{secr} 4.0 it is no longer efficient to use parallel worker processes.

}

\examples{

\dontrun{

fit0 <- list(capthist = 'captdata', model = g0~1) 
fitb <- list(capthist = 'captdata', model = g0~b)
fits <- list.secr.fit (c('fit0','fitb'))
AIC(fits)

# replacing par.derived and par.region.N:

lapply(fits, derived)
lapply(fits, region.N)

}

}

