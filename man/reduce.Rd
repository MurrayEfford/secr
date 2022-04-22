\name{reduce}
\alias{reduce}
\alias{reduce.default}

\title{ Combine Columns }

\description{
Combine columns in a matrix-like object to create a new data set using the first non-zero value.
}

\usage{
reduce (object, ...)
\method{reduce}{default}(object, columns, \dots)
}

\arguments{
  \item{object}{ object that may be coerced to a matrix }
  \item{columns}{ list in which each component is a vector of subscripts for columns to be pooled }
  \item{\dots}{ other arguments (not used currently) }
}

\details{
  
The first element of \code{columns} defines the columns of \code{object}
for the first new column, the second for the second new column etc.
This is a generic method. More useful methods exist for \code{capthist} and
\code{traps} objects.

}

\value{
A matrix with number of columns equal to \code{length(columns)}.
}

\seealso{ \code{\link{capthist}}, \code{\link{reduce.capthist}}, \code{\link{reduce.traps}} }

\examples{

## matrix with random zeros
temp <- matrix(runif(20), nc = 4)
temp[sample(20,10)] <- 0
temp

reduce(temp, list(1:2, 3:4))

}

\keyword{ manip }

