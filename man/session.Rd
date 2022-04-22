\name{session}
\alias{session}
\alias{session<-}

\title{ Session Vector }

\description{
  Extract or replace the session names of a \code{capthist} object.
}

\usage{
session(object, \dots)
session(object) <- value
}

\arguments{
  \item{object}{ object with `session' attribute e.g. \code{capthist} }
  \item{value}{ character vector or vector that may be coerced to character, one value per session }
  \item{\dots}{ other arguments (not used) }
}

\details{
  Replacement values will be coerced to character.
}

\value{
  a character vector with one value for each session in \code{capthist}
}

\note{

Like Density, \pkg{secr} uses the term `session' for a
closed-population sample. A session usually includes data from several
closely-spaced capture occasions (often consecutive days). Each `primary
session' in the `robust' design of Pollock (1982) would be treated as a
session in \pkg{secr}. \pkg{secr} also uses `session' for independent
subsets of the capture data distinguished by characteristics other than
sampling time (as above). For example, two grids trapped simultaneously
could be analysed as distinct sessions if (i) they were far enough
apart that there was negligible prospect of the same animal being caught
on both grids, and (ii) there was interest in comparing estimates from
the two grids, or fitting a common detection model.

The log likelihood for a session model is the sum of the separate
session log likelihoods. Although this assumes independence of sampling,
parameters may be shared across sessions, or session-specific parameter
values may be functions of session-level covariates. For many purposes,
`sessions' are equivalent to `groups'. For multi-session models the
detector array and mask are specified separately for each session. Group
models are therefore generally simpler to implement. On the other hand,
sessions offer more flexibility in defining and evaluating
between-session models, including trend models.

}

\seealso{ 
  \code{\link{capthist}}
}


\references{ 

Pollock, K. H. (1982) A capture-recapture design robust to unequal
probability of capture. \emph{Journal of Wildlife Management} \bold{46},
752--757.

}

\examples{
  session(captdata)
}

\keyword{ models }

