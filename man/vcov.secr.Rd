\name{vcov.secr}
\alias{vcov.secr}

\title{ Variance - Covariance Matrix of SECR Parameters }

\description{
  Variance-covariance matrix of beta or real parameters from fitted secr model.
}

\usage{
\method{vcov}{secr}(object, realnames = NULL, newdata = NULL, 
    byrow = FALSE, ...)
}

\arguments{
  \item{object}{ secr object output from the function \code{\link{secr.fit}} }
  \item{realnames}{ vector of character strings for names of `real' parameters }
  \item{newdata}{ dataframe of predictor values }
  \item{byrow}{ logical for whether to compute covariances among `real' parameters for each row of new data, or among rows for each real parameter}
  \item{\dots}{ other arguments (not used) }
}

\details{

By default, returns the matrix of variances and covariances among the
estimated model coefficients (beta parameters). 

If \code{realnames} and \code{newdata} are specified, the result is
either a matrix of variances and covariances for each `real' parameter
among the points in predictor-space given by the rows of \code{newdata}
or among real parameters for each row of \code{newdata}. Failure to
specify \code{newdata} results in a list of variances only.

}

\value{

A matrix containing the variances and covariances among beta parameters
on the respective link scales, or a list of among-parameter variance-covariance
matrices, one for each row of \code{newdata}, or a list of among-row variance-covariance
matrices, one for each `real' parameter.

}

\seealso{ \code{\link{vcov}}, \code{\link{secr.fit}}, \code{\link{print.secr}} }

\examples{
## previously fitted secr model
vcov(secrdemo.0)
}


