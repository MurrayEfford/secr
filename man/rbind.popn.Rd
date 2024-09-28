\name{rbind.popn}
\alias{rbind.popn}

\title{ Combine popn Objects }

\description{
Form a single \code{popn} object from two or more existing \code{popn} objects, or a list.
}

\usage{

\method{rbind}{popn} (..., renumber = TRUE)

}

\arguments{
  \item{\dots}{ one or more \code{popn} objects }
  \item{renumber}{ logical for whether row names in the new object should be set to the row indices }
}

\details{
  An attempt to combine objects will fail if they conflict in their \code{covariates} 
  attributes.  
  
  From \pkg{secr} 3.1 this is an S3 method and list input is not allowed.
}

\value{
An object of class \code{popn} with number of rows equal to the sum of the rows in the input objects.
}

\seealso{ \code{\link{popn}} }

\examples{

## generate and combine two subpopulations
trapobj <- make.grid()
p1 <- sim.popn(D = 3, core = trapobj)
p2 <- sim.popn(D = 2, core = trapobj)
covariates(p1) <- data.frame(size = rep("small", nrow(p1)))
covariates(p2) <- data.frame(size = rep("large", nrow(p2)))
pop <- rbind(p1,p2)

## or
pop <- do.call(rbind, list(p1,p2))

}

