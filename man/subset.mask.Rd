\name{subset.mask}
\alias{subset.mask}
\alias{split.mask}
\alias{rbind.mask}

\title{ Subset, Split or Combine Mask Objects }

\description{

  Retain selected rows of a \code{mask} object.
  
}
\usage{

\method{subset}{mask}(x, subset, ...)

\method{split}{mask}(x, f, drop = FALSE, clusters = NULL, ...)

\method{rbind}{mask}(...)
}

\arguments{
  \item{x}{ \code{mask} object }
  \item{subset}{ numeric or logical vector to select rows of mask }
  \item{f}{ factor or object that may be coerced to a factor }
  \item{drop}{ logical indicating if levels that do not occur should be dropped (if f is a factor) }
  \item{clusters}{list of traps objects, each defining a cluster (alternative to f)}
  \item{\dots}{ two or more \code{mask} objects (rbind only) }
}

\details{

The subscripts in \code{subset} may be of type integer, character or
logical as described in \code{\link{Extract}}.

The \code{split} method may use either a factor f with one value for each row or a list of clusters, each a traps object. The output mask corresponding to each cluster is the subset of the original mask points that lie within buffer of a trap within the cluster; buffer is computed as the maximum distance between a mask point in \code{x} and any detector in \code{clusters}. Sub-masks specified with \code{clusters} may overlap.

Covariates are ignored by \code{rbind.mask}.

}

\value{

  For \code{subset}, an object of class `mask' with only the requested
  subset of rows and `type' attribute set to `subset'.

  For \code{split}, a list of mask objects.

  For \code{rbind}, an object of class `mask' with all unique rows from
  the masks in \dots, and `type' attribute set to `rbind'.

}

\section{Warning}{

The spacing attribute is carried over from the input (it is not updated
automatically). In the case of very sparse masks (i.e. those with
isolated points) this may lead to an unexpected value for this
attribute. (Automatic updating requires excessive computation time
and/or memory for very large masks).
  
  }

\seealso{ \code{\link{mask}} }

\examples{
tempmask <- make.mask(make.grid())
OK <- (tempmask$x + tempmask$y) > 100
tempmask <- subset(tempmask, subset = OK)
plot(tempmask)

}


