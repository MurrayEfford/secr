\name{secr-deprecated}
\alias{secr-deprecated}
\alias{par.secr.fit}
\alias{par.region.N}
\alias{par.derived}

\title{Deprecated Functions in Package \pkg{secr}}

\description{

These functions will be removed from future versions of \pkg{secr}.

}

\usage{

par.secr.fit (arglist, ncores = 1, seed = NULL, trace = TRUE, logfile = "logfile.txt", 
    prefix = "fit.", LB = FALSE, save.intermediate = FALSE)
par.derived (secrlist, ncores = 1, ...)
par.region.N (secrlist, ncores = 1, ...)

}

\arguments{
  
  \item{arglist}{list of argument lists for \code{secr.fit} or a character
    vector naming such lists}
  
  \item{ncores}{ integer number of cores to be used for parallel processing}

  \item{seed}{integer pseudorandom number seed}

  \item{trace}{logical; if TRUE intermediate output may be logged}
  
  \item{logfile}{character name of file to log progress reports}
  
  \item{prefix}{character prefix for names of output}
  
  \item{LB}{logical; if TRUE then use load balancing}
  
  \item{save.intermediate}{logical; if TRUE then each fit is saved to an external file}

  \item{\dots}{other arguments passed to \code{derived} or \code{region.N}}

  \item{secrlist}{secrlist object}

}

\details{

Since the introduction of multi-threading in \pkg{secr} 4.0 it is no longer efficient to use parallel worker processes.

\code{\link{list.secr.fit}} replaces \code{par.secr.fit}.

Functions par.derived and par.region.N can be replaced by a simple call of lapply (see Examples in \code{\link{list.secr.fit}}).

}

\seealso{

\code{\link{secr-defunct}}, 

}

