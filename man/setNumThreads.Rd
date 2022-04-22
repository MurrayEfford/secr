\name{setNumThreads}
\alias{setNumThreads}

\title{Number of Threads}

\description{

Set or report the number of cores to be used for multi-threaded operations. A wrapper for the RcppParallel function setThreadOptions (Allaire et al. 2019).

}

\usage{

setNumThreads(ncores, ...)

}

\arguments{

  \item{ncores}{integer number of threads to use}
  \item{\dots}{other arguments passed to \code{RcppParallel::setThreadOptions}, specifically stackSize}

}

\details{

If \code{ncores} is NULL then the current value of the environment variable RCPP_PARALLEL_NUM_THREADS is used. RCPP_PARALLEL_NUM_THREADS defaults to 2 at the start of a session (assuming at least 2 logical cores available).

Calling \code{setNumThreads()} with no arguments is a handy way to check how many threads are in use.

The value of RCPP_PARALLEL_NUM_THREADS is also reset when a multi-threaded function such as \code{secr.fit} is called with a non-NULL value of the ncores argument. This value applies in later calls of \code{secr.fit} with \code{ncores = NULL} until changed.

}
\value{

The new value of the environment variable RCPP_PARALLEL_NUM_THREADS.

}

\references{

Allaire, J. J., Francois, R., Ushey, K., Vandenbrouck, G., Geelnard, M. and Intel (2019)
  RcppParallel: Parallel Programming Tools for 'Rcpp'. R package version 4.4.4.
  https://CRAN.R-project.org/package=RcppParallel.

}

\seealso{

\link{Parallel},
\code{\link{setThreadOptions}}
\code{\link{Sys.getenv}}

}

\note{

The mechanism for setting the number of threads changed between versions 4.1.0 and 4.2.0. The default number of cores is now capped at 2 to meet CRAN requirements. Setting \code{ncores = NULL} previously specified one less than the number of available cores.

}

\examples{

# determine current number of threads

setNumThreads()

\dontrun{

# set new number of threads
setNumThreads(7)

# a call to secr.fit that specifies 'ncores' also sets the 
# number of threads, as we see here

fit <- secr.fit(captdata, trace = FALSE, ncores = 8)
setNumThreads()

}

}
