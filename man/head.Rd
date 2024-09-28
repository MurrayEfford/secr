\name{head}
\alias{head.mask}
\alias{head.Dsurface}
\alias{head.traps}
\alias{head.capthist}
\alias{tail.mask}
\alias{tail.Dsurface}
\alias{tail.traps}
\alias{tail.capthist}

\title{ First or Last Part of an Object }

\description{
  Returns the first or last parts of secr objects
}

\usage{
\method{head}{mask}(x, n=6L, ...)
\method{head}{Dsurface}(x, n=6L, ...)
\method{head}{traps}(x, n=6L, ...)
\method{head}{capthist}(x, n=6L, ...)
\method{tail}{mask}(x, n=6L, ...)
\method{tail}{Dsurface}(x, n=6L, ...)
\method{tail}{traps}(x, n=6L, ...)
\method{tail}{capthist}(x, n=6L, ...)
}

\arguments{
  \item{x}{ `mask', `traps' or `capthist' object }
  \item{n}{ a single integer. If positive, size for the resulting
  object: number of elements for a vector (including lists), rows for a
  matrix or data frame or lines for a function. If negative, all but the
  n last/first number of elements of x. }
 \item{\dots}{ other arguments passed to \code{subset}}
}

\details{
  
These custom S3 methods retain the class of the target object, unlike the
default methods applied to `mask', `Dsurface', `traps' or `capthist' objects.

}

\value{
An object of the same class as x, but (usually) fewer rows.
}

\seealso{ \code{\link{head}}, \code{\link{tail}} }

\examples{
head(possummask)
}


