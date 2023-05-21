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
  \code{\link{derivedSystematic}} \cr
  \code{\link{esa}} \cr
  \code{\link{fxi.secr}} and related functions \cr
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
through the package \pkg{\link{parallel}}. This mechanism is still available in the functions listed here, but the speed gains are often small or even negative. Set \code{ncores > 1} in the function call to use multiple cores.

\tabular{llll}{
  Function \tab Unit \tab Benefit \tab Notes \cr
  \code{\link{ip.secr}} \tab replicate \tab large \tab  \cr
  \code{\link{par.secr.fit}} \tab model \tab none \tab  \cr
  \code{\link{par.derived}} \tab fitted model \tab none \tab  \cr
  \code{\link{par.region.N}} \tab fitted model \tab none \tab  \cr
}
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
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
# quad-core i7 CPU, 16 Gb RAM
# ...

## benefit from multi-threading in secr.fit

for (i in 1:8) 
     print(system.time(secr.fit(ovenCH, buffer = 400, trace = FALSE, ncores = i)))
     
#   user  system elapsed 
#  54.25    0.17   54.43 
#   user  system elapsed 
#  29.48    0.16   20.30 
#   user  system elapsed 
#  43.34    0.17   25.04 
#   user  system elapsed 
#  43.92    0.14   22.95 
#   user  system elapsed 
#  46.16    0.22   22.70 
#   user  system elapsed 
#  31.59    0.17   15.19 
#   user  system elapsed 
#  45.58    0.12   21.93 
#   user  system elapsed 
#  37.43    0.15   18.43 

## and for simulation...
 
for (i in 1:8)
     print(system.time(sim.secr(secrdemo.0, nsim = 20, tracelevel = 0, ncores = i)))
     
#   user  system elapsed 
# 160.68    0.78  161.93 
#   user  system elapsed 
# 158.48    0.56   89.42 
#   user  system elapsed 
# 158.54    0.43   65.05 
#   user  system elapsed 
# 165.00    0.39   55.06 
#   user  system elapsed 
# 171.38    0.53   48.03 
#   user  system elapsed 
# 191.03    0.47   48.14 
#   user  system elapsed 
# 184.46    0.52   43.42 
#   user  system elapsed 
# 193.07    0.56   42.34 

for (i in 1:8)
    print(system.time(LLsurface(secrdemo.0, ncores = i)))

# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   26.59    0.14   26.80 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   25.74    0.13   17.73 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   26.08    0.18   15.12 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   26.93    0.19   13.82 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   28.32    0.06   13.52 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   29.40    0.16   13.05 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   29.76    0.22   12.92 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   29.10    0.20   12.81 

}

}
