\name{print.capthist}
\alias{print.capthist}

\title{ Print Detections }

\description{
  Print method for \code{capthist} objects.
}

\usage{
\method{print}{capthist}(x, ..., condense = FALSE, sortrows = FALSE)
}

\arguments{
  \item{x}{ \code{capthist} object }
  \item{\dots}{ arguments to pass to \code{\link{print.default}} }
  \item{condense}{ logical, if true then use condensed format for 3-D data }
  \item{sortrows}{ logical, if true then sort output by animal }
}

\details{

The \code{condense} option may be used to format data from proximity
detectors in a slightly more readable form. Each row then presents the
detections of an individual in a particular trap, dropping rows (traps)
at which the particular animal was not detected.

}

\value{
Invisibly returns a dataframe (condense = TRUE) or array in the format printed.
}

\seealso{ \code{\link{print}}, \code{\link{capthist}} }

\examples{
## simulated detections of simulated default population of 5/ha
print(sim.capthist(make.grid(nx=5,ny=3)))
}


