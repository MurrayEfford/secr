\name{secr-deprecated}
\alias{secr-deprecated}
\alias{par.secr.fit}
\alias{par.region.N}
\alias{par.derived}
\alias{fxi.contour}
\alias{fxi.mode}
\alias{fx.total}
\alias{esa.plot}

\title{Deprecated Functions in Package \pkg{secr}}

\description{

These functions will be removed from future versions of \pkg{secr} or have been renamed.

par.secr.fit (arglist, ncores = 1, seed = NULL, trace = TRUE, logfile = "logfile.txt", 
    prefix = "fit.", LB = FALSE, save.intermediate = FALSE)

par.derived (secrlist, ncores = 1, ...)

par.region.N (secrlist, ncores = 1, ...)

fxi.contour (object, i = 1, sessnum = 1, border = 100, nx = 64, levels = NULL, 
    p = seq(0.1,0.9,0.1), plt = TRUE, add = FALSE, fitmode = FALSE, plotmode = 
    FALSE, fill = NULL, output = c('list','sf','SPDF'), ncores = NULL, ...) 

fxi.mode (object, i = 1, sessnum = 1, start = NULL, ncores = NULL, ...)

fx.total (object, sessnum = 1, mask = NULL, ncores = NULL, ...)

esa.plot (object, max.buffer = NULL, spacing = NULL, max.mask = NULL, detectfn,
      detectpar, noccasions, binomN = NULL, thin = 0.1, poly = NULL,
      poly.habitat = TRUE, session = 1, plt = TRUE, type = c('density', 
      'esa','meanpdot', 'CVpdot'), n = 1, add = FALSE, overlay = TRUE, 
      conditional = FALSE, ...)

}

\details{

Since the introduction of multi-threading in \pkg{secr} 4.0 it is no longer efficient to use parallel worker processes.

\code{\link{list.secr.fit}} replaces \code{par.secr.fit}.

Functions par.derived and par.region.N can be replaced by a simple call of lapply (see Examples in \code{\link{list.secr.fit}}).

Some functions have been renamed to avoid the ambiguous ".":

\tabular{lll}{
Old \tab New \cr
fxi.secr \tab \code{\link{fxi}} (S3 generic with method for secr objects) \cr
fxi.contour \tab \code{\link{fxiContour}} \cr
fxi.mode \tab \code{\link{fxiMode}} \cr
fx.total \tab \code{\link{fxTotal}} \cr
esa.plot \tab esaPlot \cr
}

}

\seealso{

\code{\link{secr-defunct}}, \code{\link{list.secr.fit}}

}

