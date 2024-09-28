\name{clone}
\alias{clone}
\alias{clone.default}
\alias{clone.popn}
\alias{clone.capthist}

\title{
Replicate Rows
}

\description{
  
  Clone rows of an object a constant or random number of times

}

\usage{
  \method{clone}{default}(object, type, ...)
  \method{clone}{popn}(object, type, ...)
  \method{clone}{capthist}(object, type, ...)
}

\arguments{
  \item{object}{any object}
  \item{type}{character `constant', `poisson',  `truncatedpoisson' or `nbinom'}
  \item{\dots}{other arguments for distribution function }
}
\details{
  
The \dots argument specifies the number of times each row should be
repeated. For random distributions (Poisson or negative binomial) \dots
provides the required parameter values: \code{lambda} for Poisson,
\code{size, prob} or \code{size, mu} for negative binomial.

One application is to derive a population of cues from a popn object,
where each animal in the original popn generates a number of cues from
the same point.

Cloning a capthist object replicates whole detection
histories. Individual covariates and detection-specific attributes
(e.g., signal strength or xy location in polygon) are also
replicated. Cloned data from single-catch traps will cause verify() to
fail, but a model may still be fitted in \code{secr.fit} by overriding
the check with \code{verify = FALSE}.

}

\value{

Object of same class as \code{object} but with varying number of
rows. For \code{clone.popn} and \code{capthist} an attribute `freq' is
set, a vector of length equal to the original number of rows giving the
number of repeats (including zeros).

If \code{popn} or \code{capthist} is a multi-session object the returned value will be
a multi-session object of the same length.

}

\seealso{
\code{\link{sim.popn}}
}

\examples{

## population of animals at 1 / hectare generates random
## Poisson number of cues, lambda = 5
mics4 <- make.grid( nx = 2, ny = 2, spacing = 44, detector = "signal")
pop <- sim.popn (D = 1, core = mics4, buffer = 300, nsessions = 6)
pop <- clone (pop, "poisson", 5)
attr(pop[[1]],"freq")

clone(captdata, "poisson", 3)

# To avoid losing any individuals use zero-truncated Poisson
# First find lambda of truncated Poisson with given mean
getlambda <- function (target) {
    fn <- function(x) x / (1-exp(-x)) - target
    uniroot(interval = c(1e-8, target), f = fn)$root
}
clone(captdata, "truncatedpoisson", getlambda(3))

}

