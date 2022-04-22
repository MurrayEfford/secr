\name{housemouse}
\docType{data}

\alias{housemouse}
\alias{Coulombe}

\title{ House mouse live trapping data  }

\description{
Data of H. N. Coulombe from live trapping of feral house mice (\emph{Mus
 musculus}) in a salt marsh, California, USA.
}

\usage{housemouse}

\details{

H. N. Coulombe conducted a live-trapping study on an outbreak of feral house
mice in a salt marsh in mid-December 1962 at Ballana Creek, Los Angeles
County, California. A square 10 x 10 grid was used with 100 Sherman
traps spaced 3 m apart. Trapping was done twice daily, morning and
evening, for 5 days.

The dataset was described by Otis et al. (1978) and distributed with
their CAPTURE software (now available from
\url{https://www.mbr-pwrc.usgs.gov/software.html}). Otis et
al. (1978 p. 62, 68) cite Coulombe's unpublished 1965 master's thesis
from the University of California, Los Angeles, California.

The data are provided as a single-session \code{capthist} object. There
are two individual covariates: sex (factor levels `f', `m') and age
class (factor levels `j', `sa', `a'). The sex of two animals is not
available (NA); it is necessary to drop these records for analyses
using `sex' unless missing values are specifically allowed, as in \code{\link{hcov}}.

The datasets were originally in the CAPTURE `xy complete' format which
for each detection gives the `column' and `row' numbers of the trap
(e.g. ` 9 5' for a capture in the trap at position (x=9, y=5) on the
grid). Trap identifiers have been recoded as strings with no spaces by
inserting zeros (e.g. `0905' in this example).

Sherman traps are designed to capture one animal at a time, but the data
include 30 double captures and one occasion when there were 4
individuals in a trap at one time. The true detector type therefore
falls between `single' and `multi'. Detector type is set to `multi' in
the distributed data objects.

Otis et al. (1978) report various analyses including a closure test on
the full data, and model selection and density estimation on data from
the mornings only. 

}

\source{

File `examples' distributed with program CAPTURE.

}

\references{

Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978)
  Statistical inference from capture data on closed animal
  populations. \emph{Wildlife Monographs} \bold{62}, 1--135.

}

\examples{

plot(housemouse, title = paste("Coulombe (1965), Mus musculus,",
    "California salt marsh"), border = 5, rad = 0.5,
    gridlines = FALSE)

morning <- subset(housemouse, occ = c(1,3,5,7,9))
summary(morning)

## drop 2 unknown-sex mice
known.sex <- subset(housemouse, !is.na(covariates(housemouse)$sex))

## reveal multiple captures
table(trap(housemouse), occasion(housemouse))

\dontrun{

## assess need to distinguish morning and afternoon samples
housemouse.0 <- secr.fit (housemouse, buffer = 20)
housemouse.ampm <- secr.fit (housemouse, model = g0~tcov, buffer = 20,
    timecov = c(0,1,0,1,0,1,0,1,0,1))
AIC(housemouse.0, housemouse.ampm)

}

}

\keyword{datasets}
