\name{rbind.traps}
\alias{rbind.traps}

\title{ Combine traps Objects }

\description{
Form a single \code{traps} object from two or more existing \code{traps} objects.
}

\usage{
\method{rbind}{traps}(..., renumber = TRUE, addusage, checkdetector = TRUE, suffix = TRUE)
}
\arguments{
  \item{\dots}{ one or more \code{traps} objects }
  \item{renumber}{ logical for whether row names in the new object
    should be set to the row indices }
  
  \item{addusage}{ integer vector; if specified and the inputs lack
    usage attributes then a binary usage attribute will be generated
    with the given number of occasions for each input }
    \item{checkdetector}{logical; if TRUE then variation in the detector attribute triggers a warning }
    \item{suffix}{logical; if TRUE then suffix to the row names indicates source }

}

\details{
  
  An attempt to combine objects will fail if they conflict in their
  \code{covariates} attributes. Differences in the \code{usage}
  attribute are handled as follows. If \code{usage} is missing for all
  inputs and \code{addusage = TRUE} is specified then usage codes are
  generated automatically (positive for the specified number of
  occasions). If \code{usage} is specified for one input but not
  other(s), the missing values are constructed assuming all detectors
  were operated for the maximum number of occasions in any input. If
  inputs differ in the number of `usage' columns (occasions), the
  smaller matrices are padded with `zero' columns to the maximum number
  of columns in any input.
  
  \dots may be a single multi-session traps object (from 2.10.0).
  
  By default (and always prior to 3.1.1) row names include a suffix 
  (e.g., ".1", or ".2") to indicate the original object (first, second etc.). 
  A suffix is added automatically to all names if any name is duplicated, and a warning is generated.

}

\value{

An object of class \code{traps} with number of rows equal to the sum of the rows in the input objects.

}

\seealso{ \code{\link{traps}}, \code{\link{subset.traps}} }

\examples{

## nested hollow grids
hollow1 <- make.grid(nx = 8, ny = 8, hollow = TRUE)
hollow2 <- shift(make.grid(nx = 6, ny = 6, hollow = TRUE), 
    c(20, 20))
nested <- rbind (hollow1, hollow2)
plot(nested, gridlines = FALSE, label = TRUE)
}

\keyword{ manip }
