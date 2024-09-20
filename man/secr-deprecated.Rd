\name{secr-deprecated}
\alias{secr-deprecated}
\alias{par.secr.fit}
\alias{par.region.N}
\alias{par.derived}
\alias{fxi.contour}
\alias{fxi.mode}
\alias{fxi.secr}
\alias{fx.total}

\title{Deprecated Functions in Package \pkg{secr}}

\description{

These functions will be removed from future versions of \pkg{secr} or have been renamed.

}

\usage{

par.secr.fit (arglist, ncores = 1, seed = NULL, trace = TRUE, logfile = "logfile.txt", 
    prefix = "fit.", LB = FALSE, save.intermediate = FALSE)
par.derived (secrlist, ncores = 1, ...)
par.region.N (secrlist, ncores = 1, ...)

fxi.contour (object, i = 1, sessnum = 1, border = 100, nx = 64, levels = NULL, 
    p = seq(0.1,0.9,0.1), plt = TRUE, add = FALSE, fitmode = FALSE, plotmode = 
    FALSE, fill = NULL, output = c('list','sf','SPDF'), ncores = NULL, ...) 
fxi.mode (object, i = 1, sessnum = 1, start = NULL, ncores = NULL, ...)
fx.total (object, sessnum = 1, mask = NULL, ncores = NULL, ...)
fxi.secr (object, i = NULL, sessnum = 1, X = NULL, ncores = NULL) 

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

Some functions will be renamed to avoid the ambiguous "." --

\tabular{lll}{
Old \tab New \cr
fxi.contour \tab fxiContour \cr
fxi.mode \tab fxiMode \cr
fxi.secr \tab fxi \cr
fx.total \tab fxTotal \cr
}

}

\seealso{

\code{\link{secr-defunct}}, 

}

