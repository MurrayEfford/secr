\name{covariates}
\alias{covariates}
\alias{covariates<-}

\title{ Covariates Attribute }

\description{
  Extract or replace covariates
}

\usage{
covariates(object, ...)
covariates(object) <- value
}

\arguments{
  \item{object}{ an object of class \code{traps}, \code{popn}, \code{capthist}, or \code{mask} }
  \item{value}{ a dataframe of covariates }
  \item{\dots}{ other arguments (not used) }
}

\details{
  For replacement, the number of rows of \code{value} must match exactly the number of rows in \code{object}.
}

\value{
  covariates(object) returns the dataframe of covariates associated with
  \code{object}. \code{covariates(object)} may be NULL.

  Individual covariates are stored in the `covariates' attribute of a
  \code{capthist} object.

  Covariates used for modelling density are stored in the `covariates'
  attribute of a \code{mask} object.
  
  Detector covariates may vary between sampling occasions. In this case,
   columns in the detector covariates data.frame are associated with
   particular times; the matching is controlled by the
   \code{\link{timevaryingcov}} attribute.
}

\seealso{

  \code{\link{timevaryingcov}}
  
}

\examples{
## detector covariates
temptrap <- make.grid(nx = 6, ny = 8)
covariates (temptrap) <- data.frame(halfnhalf = 
    factor(rep(c("left","right"),c(24,24))) )
summary(covariates(temptrap))
}

\keyword{ manip }

