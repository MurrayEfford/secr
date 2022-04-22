\name{signalmatrix}
\alias{signalmatrix}

\title{
Reformat Signal Data
}

\description{
Produce sound x microphone matrix, possibly with sound covariates as
extra columns.
}

\usage{
signalmatrix(object, noise = FALSE, recodezero = FALSE,
    prefix = "Ch", signalcovariates = NULL, names = NULL)

}

\arguments{
  \item{object}{object inheriting from secr class `capthist'}
  \item{noise}{logical; if TRUE, noise is extracted instead of signal}
  \item{recodezero}{logical; if TRUE zero signals are set to NA}
  \item{prefix}{character value used to form channel names}
  \item{signalcovariates}{character vector of covariate names from signalframe to
    add as columns}
  \item{names}{character vector of column names}
}

\details{
  
This function extracts signal or noise data from a capthist object, where
it is stored in the `signalframe' attribute, as a more natural sound x microphone
table. There is no equivalent replacement function.

The \code{signalcovariates} argument may be used to specify additional
columns of the signal frame to collapse and add as columns to the right
of the actual signal data. Ordinarily there will be multiple rows in
signalframe for each row in the output; the covariate value is taken
from the first matching row.

If \code{names} is not provided, column names are constructed from the
detector names. If the length of \code{names} is less than the number of
columns, simple numerical names are constructed.

}

\value{

A dataframe with dim = c(n,K+j) where n is the
number of separate sounds, K is the number of microphones, and j is the
number of covariates (by default j = 0).

}

\seealso{
\code{\link{ovensong}}
}

\examples{

## use 'secr' ovenbird data
signalmatrix(signalCH)

}
\keyword{ manip }
