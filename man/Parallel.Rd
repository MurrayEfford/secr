\name{Parallel}
\alias{Parallel}
\alias{ncores}
\alias{Multi-core processing}

\title{ Multi-core Processing }

\description{

From version 4.0 \pkg{secr} uses multi-threading in C++ (package \pkg{RcppParallel}, Allaire et al. 2021) to speed likelihood evaluation and hence model fitting in \code{secr.fit}. Detection histories are distributed over threads. Setting \code{ncores = NULL} in functions with multi-threading uses the existing value from the environment variable RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

These functions use multi-threading and call \code{setNumThreads} internally:
\tabular{llll}{
  \code{\link{autoini}} \cr
  \code{\link{confint.secr}} \cr
  \code{\link{derived.secr}} \cr
  \code{\link{esa}} \cr
  \code{\link{fxi}} and related functions \cr
  \code{\link{pdot}} \cr
  \code{\link{region.N}} \cr
  \code{\link{score.test}} \cr
  \code{\link{secr.fit}} \cr
}


These functions use multi-threading without calling \code{setNumThreads}:
\tabular{llll}{
  \code{\link{LLsurface.secr}} \cr
  \code{\link{mask.check}} \cr
  \code{\link{expected.n}} \cr
  \code{\link{secr.test}} \cr
  \code{\link{sim.secr}} \cr
}

Other functions may use multithreading indirectly through a call to one of these functions. Examples are \code{\link{suggest.buffer}} (\code{autoini}), \code{\link{esa.plot}} (\code{pdot}), and \code{\link{bias.D}} (\code{pdot}).

NOTE: The mechanism for setting the number of threads changed between versions 4.1.0 and 4.2.0. The default number of cores is now capped at 2 to meet CRAN requirements. Setting \code{ncores = NULL} previously specified one less than the maximum number of cores.

Earlier versions of \pkg{secr} made more limited use of multiple cores (CPUs)
through the package \pkg{\link{parallel}}. The functions \code{\link{par.secr.fit}}, \code{\link{par.derived}}, and   \code{\link{par.region.N}} are now deprecated because they were too slow. \code{\link{list.secr.fit}} replaces \code{par.secr.fit}

`Unit' refers to the unit of work sent to each worker process. As a guide, a `large' benefit means >60\% reduction in process time with
  4 CPUs.

\pkg{parallel} offers several different mechanisms, bringing together
the functionality of \pkg{multicore} and \pkg{snow}. The mechanism used
by \pkg{secr} is the simplest available, and is expected to work across all
operating systems. Technically, it relies on Rscript and communication
between the master and worker processes is \emph{via} sockets. As stated
in the \pkg{parallel} documentation "Users of Windows and Mac OS X may
expect pop-up dialog boxes from the firewall asking if an R process
should accept incoming connections". You may possibly get warnings from R 
about closing unused connections. These can safely be ignored.

Use \code{parallel::detectCores()} to get
an idea of how many cores are available on your machine; this may (in
Windows) include virtual cores over and above the number of physical
cores. See RShowDoc("parallel", package = "parallel") in core R for
explanation.

In \code{secr.fit} the output component `proctime' misrepresents the
elapsed processing time when multiple cores are used.

}

\references{

Allaire, J. J., Francois, R., Ushey, K., Vandenbrouck, G., Geelnard, M. and Intel (2021)
  RcppParallel: Parallel Programming Tools for 'Rcpp'. R package version 5.1.2.
  https://CRAN.R-project.org/package=RcppParallel.

}

\section{Warning}{

It appears that multicore operations in \pkg{secr} using \pkg{parallel} may fail if the packages \pkg{snow} and \pkg{snowfall} are also loaded. The error message is obscure:

``Error in UseMethod("sendData") : 
  no applicable method for 'sendData' applied to an object of class "SOCK0node"''

}

\examples{

\dontrun{

sessionInfo()
# R version 4.3.0 (2023-04-21 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22621)
# on Dell-XPS 8950 Intel i7-12700K
# ...
# see stackoverflow suggestion for microbenchmark list argument
# https://stackoverflow.com/questions/32950881/how-to-use-list-argument-in-microbenchmark

library(microbenchmark)
options(digits = 4)

## benefit from multi-threading in secr.fit

jobs <- lapply(seq(2,8,2), function(nc) 
    bquote(suppressWarnings(secr.fit(captdata, trace = FALSE, ncores = .(nc)))))
microbenchmark(list = jobs, times = 10, unit = "seconds")
# [edited output]
# Unit: seconds
# expr     min      lq   mean median     uq    max neval
# ncores = 2 1.75880 2.27978 2.6680 2.7450 3.0960 3.4334    10
# ncores = 4 1.13549 1.16280 1.6120 1.4431 2.0041 2.4018    10
# ncores = 6 0.88003 0.98215 1.2333 1.1387 1.5175 1.6966    10
# ncores = 8 0.78338 0.90033 1.5001 1.0406 1.2319 4.0669    10

## sometimes (surprising) lack of benefit with ncores>2

msk <- make.mask(traps(ovenCH[[1]]), buffer = 300, nx = 25)
jobs <- lapply(c(1,2,4,8), function(nc) 
    bquote(secr.fit(ovenCH, trace = FALSE, ncores = .(nc), mask = msk)))
microbenchmark(list = jobs, times = 10, unit = "seconds")
# [edited output]
# Unit: seconds
# expr     min      lq   mean median     uq    max neval
# ncores = 1 12.5010 13.4951 15.674 15.304 16.373 21.723    10
# ncores = 2 10.0363 11.8634 14.407 13.726 16.782 22.966    10
# ncores = 4  8.6335 10.3422 13.085 12.449 15.729 17.914    10
# ncores = 8  8.5286  9.9008 10.751 10.736 10.796 14.885    10

## and for simulation...

jobs <- lapply(seq(2,8,2), function(nc)
    bquote(sim.secr(secrdemo.0, nsim = 20, tracelevel = 0,  ncores = .(nc))))
microbenchmark(list = jobs, times = 10, unit = "seconds")
# [edited output]
# Unit: seconds
# expr    min     lq   mean median     uq     max neval
# ncores = 2 48.610 49.932 59.032 52.485 54.730 119.905    10
# ncores = 4 29.480 29.996 30.524 30.471 31.418  31.612    10
# ncores = 6 22.583 23.594 24.148 24.354 24.644  25.388    10
# ncores = 8 19.924 20.651 25.581 21.002 21.696  51.920    10

## and log-likelihood surface

jobs <- lapply(seq(2,8,2), function(nc) 
    bquote(suppressMessages(LLsurface(secrdemo.0,  ncores = .(nc)))))
microbenchmark(list = jobs, times = 10, unit = "seconds")
# [edited output]
# Unit: seconds
# expr    min     lq   mean median     uq    max neval
# ncores = 2 20.941 21.098 21.290 21.349 21.471 21.619    10
# ncores = 4 14.982 15.125 15.303 15.263 15.449 15.689    10
# ncores = 6 13.994 14.299 14.529 14.342 14.458 16.515    10
# ncores = 8 13.597 13.805 13.955 13.921 14.128 14.353    10

}

}
