\name{subset.traps}
\alias{subset.traps}
\alias{split.traps}

\title{ Subset traps Object }

\description{
  Retain selected rows of a traps object.
}
\usage{
\method{subset}{traps}(x, subset = NULL, occasions = NULL, \dots)
\method{split}{traps}(x, f, drop = FALSE, prefix = "S", byoccasion = FALSE, \dots)
}

\arguments{
  \item{x}{ \code{traps} object }
  \item{subset}{ vector to subscript the rows of \code{x} }
  \item{occasions}{ vector to subscript columns in \code{usage(x)}}
  \item{\dots}{ arguments passed to other functions or to optional subset function (subset.traps)  }
  \item{f}{ factor or object that may be coerced to a factor }
  \item{drop}{ logical indicating if levels that do not occur should be dropped (if f is a factor) }
  \item{prefix}{ a character prefix to be used for component names when
    values of f are numeric }
  \item{byoccasion}{logical ; if TRUE then f is used to split occasions}
}
\details{

The subscripts in \code{subset} may be of type integer, character or
logical as described in \code{\link{Extract}}. By default, all rows are retained.

In the case of `polygon' and `transect' detectors, subsetting is done at
the level of whole polygons or transects. \code{subset} should therefore
have the same length as \code{levels(polyID(x))} or
\code{levels(transectID(x))}.

\code{split} generates a list in which each component is a \code{traps}
object. Each component corresponds to a level of \code{f}.  The argument
\code{x} of \code{split} cannot be a list (i.e. \code{x} must be a
single-session traps object). 

If the levels of \code{f} are numeric, from version 2.10.3 a leading zero is inserted in the names of the output list to maintain the sort order.

}
\value{

An object of class \code{traps} with only the requested subset of rows.
Subsetting is applied to \code{usage} and \code{covariates} attributes
if these are present.

Splitting with \code{byoccasion = TRUE} produces a list of traps
objects, each with usage codes for a subset of occasions. Traps not used
on any occasion within a session are automatically dropped from that
session.

}

\seealso{ 

\code{\link{traps}}, 
\code{\link{rbind.traps}} 

}

\section{Warning}{

\code{split.traps} does not work for mark--resight data.

}

\examples{
## odd-numbered traps only, using modulo operator
temptrap <- make.grid(nx = 7, ny = 7)
t2 <- subset(temptrap, as.logical(1:nrow(temptrap) \%\% 2))
plot(t2)

## this works also for even number of rows, but must change 'outer' call
temptrap <- make.grid(nx = 8, ny = 8)
t3 <- subset(temptrap, !as.logical(outer(1:8,1:8,'+')\%\%2))
plot(t3)
}

\keyword{ manip }

