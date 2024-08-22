\name{deermouse}
\docType{data}

\alias{deermouse}
\alias{deermouse.ESG}
\alias{deermouse.WSG}

\title{ Deermouse Live-trapping Datasets }

\description{

Data of V. H. Reid from live trapping of deermice (\emph{Peromyscus
  maniculatus}) at two sites in Colorado, USA.

}

\usage{

deermouse.ESG
deermouse.WSG

}

\details{

Two datasets of V. H. Reid were described by Otis et al. (1978) and
distributed with their CAPTURE software (now available from
\url{https://eesc.usgs.gov/mbr/software/capture.shtml}). They have been used in
several other papers on closed population methods (e.g., Huggins 1991,
Stanley and Richards 2005). This description is based on pages 32 and
87--93 of Otis et al. (1978).

Both datasets are from studies in Rio Blanco County, Colorado, in the
summer of 1975. Trapping was for 6 consecutive nights. Traps were
arranged in a 9 x 11 grid and spaced 50 feet (15.2 m) apart.

The first dataset was described by Otis et al. (1978: 32) as from `a
drainage bottom of sagebrush, gambel oak, and serviceberry with pinyon
pine and juniper on the uplands'. By matching with the `examples' file
of CAPTURE this was from East Stuart Gulch (ESG).

The second dataset (Otis et al. 1978: 87) was from Wet Swizer Creek or
Gulch (WSG) in August 1975. No specific vegetation description is given
for this site, but it is stated that Sherman traps were used and
trapping was done twice daily.

Two minor inconsistencies should be noted. Although Otis et al. (1978)
said they used data from morning trap clearances, the capture histories
in `examples' from CAPTURE include a `pm' tag on each record. We
assume the error was in the text description, as their numerical results
can be reproduced from the data file. Huggins (1991) reproduced the East
Stuart Gulch dataset (omitting spatial data that were not relevant to his
method), but omitted two capture histories.

The data are provided as two single-session \code{capthist} objects
`deermouse.ESG' and `deermouse.WSG'. Each has a dataframe of individual
covariates, but the fields differ between the two study areas. The
individual covariates of deermouse.ESG are sex (factor levels `f', `m'),
age class (factor levels `y', `sa', `a') and body weight in grams. The
individual covariates of deermouse.WSG are sex (factor levels `f',`m')
and age class (factor levels `j', `y', `sa', `a') (no data on body
weight). The aging criteria used by Reid are not recorded.

The datasets were originally in the CAPTURE `xy complete' format which
for each detection gives the `column' and `row' numbers of the trap
(e.g. ` 9 5' for a capture in the trap at position (x=9, y=5) on the
grid). Trap identifiers have been recoded as strings with no spaces by
inserting zeros (e.g. `905' in this example).

Sherman traps are designed to capture one animal at a time, but the data
include double captures (1 at ESG and 8 at WSG -- see Examples). The true
detector type therefore falls between `single' and `multi'. Detector
type is set to `multi' in the distributed data objects.

\tabular{ll}{
  Object \tab Description \cr
 deermouse.ESG \tab capthist object, East Stuart Gulch \cr
 deermouse.WSG \tab capthist object, Wet Swizer Gulch \cr
}

}

\source{

File `examples' distributed with program CAPTURE.

}

\references{

Huggins, R. M. (1991) Some practical aspects of a conditional likelihood
  approach to capture experiments. \emph{Biometrics} \bold{47}, 725--732.

Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978)
  Statistical inference from capture data on closed animal
  populations. \emph{Wildlife Monographs} \bold{62}, 1--135.

Stanley, T. R. and Richards, J. D. (2005) A program for testing
  capture--recapture data for closure. \emph{Wildlife Society Bulletin}
  \bold{33}, 782--785.

}

\seealso{

\code{\link{closure.test}}

}

\examples{

par(mfrow = c(1,2), mar = c(1,1,4,1))
plot(deermouse.ESG, title = "Peromyscus data from East Stuart Gulch",
    border = 10, gridlines = FALSE, tracks = TRUE)
plot(deermouse.WSG, title = "Peromyscus data from Wet Swizer Gulch",
    border = 10, gridlines = FALSE, tracks = TRUE)

closure.test(deermouse.ESG, SB = TRUE)

## reveal multiple captures
table(trap(deermouse.ESG), occasion(deermouse.ESG))
table(trap(deermouse.WSG), occasion(deermouse.WSG))

}

\keyword{datasets}
