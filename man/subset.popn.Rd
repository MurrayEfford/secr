\name{subset.popn}
\alias{subset.popn}

\title{ Subset popn Object }

\description{
  Retain selected rows of a popn object.
}

\usage{
\method{subset}{popn}(x, subset = NULL, sessions = NULL, poly = NULL,
    poly.habitat = TRUE, keep.poly = TRUE, renumber = FALSE, \dots)
}

\arguments{
  \item{x}{ \code{popn} object }
  \item{subset}{ vector to subscript the rows of \code{x} }
  \item{sessions}{ vector to subscript sessions if \code{x} is a multi-session
    population }
  \item{poly}{ bounding polygon (see Details)}
  \item{poly.habitat}{ logical for whether poly represents habitat or
  its inverse (non-habitat)}
  \item{keep.poly}{ logical; if TRUE any bounding polygon is saved as the
  attribute `polygon' }
  \item{renumber}{ logical for whether to renumber rows in output }
  \item{\dots}{ arguments passed to other functions }
}
\details{

The subscripts in \code{subset} may be of type integer, character or
logical as described in \code{\link{Extract}}. By default, all rows are
retained.

In the case of a multi-session popn object (a list of populations),
\code{subset} may be a list with one component for the subscripts in
each new session.

If \code{poly} is specified, points outside \code{poly} are
dropped. \code{poly} may be one of the types descrbed in 
\code{\link{boundarytoSF}}.

}

\value{

An object of class \code{popn} with only the requested subset of rows.
Subsetting is applied to the \code{covariates} attribute if this is
present. Attributes `Ndist' and `model2D' are set to NULL.

}

\seealso{ \code{\link{popn}} }

\examples{

temppop <- sim.popn (D = 10, expand.grid(x = c(0,100), y =
    c(0,100)), buffer = 50)
## 50\% binomial sample of simulated population
temppops <- subset(temppop, runif(nrow(temppop)) < 0.5)
plot(temppop)
plot(temppops, add = TRUE, pch = 16)

}

\keyword{ manip }

