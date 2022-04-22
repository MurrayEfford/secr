\name{par.secr.fit}
\alias{par.secr.fit}
\alias{par.derived}
\alias{par.region.N}

\title{Fit Multiple SECR Models}

\description{
  
  These functions are wrappers for \code{\link{secr.fit}}, \code{\link{region.N}} and
  \code{\link{derived}} that allow multiple models to be fitted
  simultaneously on multi-core machines (but see Warning).

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

  From version 4.0, an \code{arglist} may specify \code{ncores} for a
  particular secr fit (previously ncores was ignored in \code{arglist}). 

  \code{trace} overrides any settings in \code{arglist}. Reporting of
  intermediate results is unreliable on Windows when \code{ncores > 1}.

  It is convenient to provide the names of the capthist and mask arguments
  in each component of arglist as character values (i.e. in quotes); objects thus
  named are exported from the workspace to each worker process (see Examples).

  Setting \code{LB = TRUE} when \code{ncores} > 1 causes the function to call \code{\link{clusterApplyLB}} instead of \code{\link{clusterApply}}. Load balancing in \code{clusterApplyLB} is likely to result in faster completion than the default if fits differ in their in execution time and \code{ncores < length(arglist)}, but this cannot be guaranteed owing to the additional communication required with the worker processes. Results with \code{LB = TRUE} for a given \code{seed} may not be reproducible.
  
  \code{save.intermediate} causes each fit to be saved to a file with extension .RData.

}

\value{

For \code{par.secr.fit} - \code{\link{secrlist}} of model fits (see
\code{\link{secr.fit}}). Names are created by prefixing \code{prefix} to the
names of \code{argslist}. If \code{trace} is TRUE then the total
execution time and finish time are displayed.

For \code{par.derived} - a list of dataframes output from
\code{\link{derived}}, applied to each model in turn.

For \code{par.region.N} - a list of dataframes output from
\code{\link{region.N}}, applied to each model in turn.

}

\section{Warning}{

With the introduction of multi-threading in \pkg{secr} 4.0, \code{par.secr.fit} has lost its speed advantage.

}

\seealso{
  \code{\link{secr.fit}},
  \code{\link{region.N}},
  \code{\link{derived}},
  \code{\link{AIC.secr}},
  \code{\link{predict.secr}},
  \link{Parallel}
}

\examples{

\dontrun{

fit0 <- list(capthist = 'captdata', model = g0~1) 
fitb <- list(capthist = 'captdata', model = g0~b)
fits <- par.secr.fit (c('fit0','fitb'))
AIC(fits)

par.derived(fits, se.esa = FALSE)

par.region.N(fits)

}

}
\keyword{ model }

