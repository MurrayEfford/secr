\name{skink}
\docType{data}

\alias{skink}
\alias{infraCH}
\alias{lineoCH}
\alias{LStraps}

\title{ Skink Pitfall Data }

\description{

Data from a study of skinks (\emph{Oligosoma infrapunctatum} and
\emph{O. lineoocellatum}) in New Zealand.

}

\usage{

infraCH
lineoCH
LStraps

}

\details{

Lizards were studied over several years on a steep bracken-covered
hillside on Lake Station in the Upper Buller Valley, South Island, New
Zealand. Pitfall traps (sunken cans baited with a morsel of fruit in
sugar syrup) were set in two large grids, each 11 x 21 traps nominally 5
meters apart, surveyed by tape and compass (locations determined later
with precision surveying equipment - see Examples). Three diurnal lizard
species were trapped: \emph{Oligosoma infrapunctatum},
\emph{O. lineoocellatum} and \emph{O. polychroma} (Scincidae). The
smallest species \emph{O. polychroma} was seldom caught and these data
are not included. The two other species are almost equal in average size
(about 160 mm total length); they are long-lived and probably mature in
their second or third year. The study aimed to examine their habitat use
and competitive interactions.

Traps were set for 12 3-day sessions over 1995--1996, but some sessions
yielded very few captures because skinks were inactive, and some sessions were
incomplete for logistical reasons. The data are from sessions 6 and 7 in
late spring (17--20 October 1995 and 14--17 November 1995). Traps were
cleared daily; the few skinks present when traps were closed on the
morning of the fourth day are treated as Day 3 captures. Individuals
were marked uniquely by clipping one toe on each foot. Natural toe loss
caused some problems with long-term identification; captures were
dropped from the dataset when identity was uncertain. Released animals
were occasionally recaptured in a different trap on the same day; these
records were also discarded.

The data are provided as two two-session \code{capthist} objects
`infraCH' and `lineoCH'. Also included is `LStraps', the \code{traps}
object with the coordinates and covariates of the trap sites (these data
are also embedded in each of the \code{capthist} objects). Pitfall traps
are multi-catch traps so \code{detector(LStraps)} = `multi'.  

Habitat data for each trap site are included as a dataframe of trap
covariates in \code{LStraps}. Ground cover and vegetation were recorded
for a 1-m radius plot at each trap site. The dataframe also gives the
total number of captures of each species by site on 31 days between
April 1995 and March 1996, and the maximum potential annual solar
radiation calculated from slope and aspect (Frank and Lee 1966). Each
site was assigned to a habitat class by fuzzy clustering (Kaufman
and Rousseauw 1990; package \pkg{cluster}) of a distance matrix using
the ground cover, vegetation and solar radiation variables. Sites in
class 1 were open with bare ground or low-canopy vegetation including
the heath-like \emph{Leucopogon fraseri} and grasses; sites in class 2
had more-closed vegetation, lacking \emph{Leucopogon fraseri} and with a
higher canopy that often included \emph{Coriaria arborea}. Site
variables are listed with definitions in the attribute
\code{habitat.variables} of \code{LStraps} (see Examples).

\if{html}{
  \tabular{ll}{
  Object \tab Description \cr
  infraCH \tab multi-session capthist object O. infrapunctatum \cr
  lineoCH \tab multi-session capthist object O. lineoocellatum \cr
  LStraps \tab traps object -- Lake Station grids \cr
}
}
\if{latex}{
\tabular{ll}{
  Object \tab Description \cr
  infraCH \tab multi-session capthist object \emph{O. infrapunctatum} \cr
  lineoCH \tab multi-session capthist object \emph{O. lineoocellatum} \cr
  LStraps \tab traps object -- Lake Station grids \cr
}
}

}

\source{

M. G. Efford, B. W. Thomas and N. J. Spencer unpublished data.

}

\seealso{

  \code{\link{capthist}}, \code{\link{covariates}}

}

\references{

Efford, M. G., Spencer, N. J., Thomas, B. W., Mason, R. F. and Williams,
P. In prep. Distribution of sympatric skink species in relation to habitat.

Frank, E. C. and Lee , R. (1966) Potential solar beam irradiation on
slopes. \emph{United States Forest Service Research Paper} RM-118.

Kaufman, L. and Rousseauw, P. J. (1990) \emph{Finding groups in data: an
introduction to cluster analysis.} John Wiley & Sons, New York.

Spencer, N. J., Thomas, B. W., Mason, R. F. and Dugdale, J. S. (1998)
Diet and life history variation in the sympatric lizards \emph{Oligosoma
nigriplantare polychroma} and \emph{Oligosoma lineoocellatum}. New
Zealand Journal of Zoology 25: 457--463.

}

\examples{
summary (infraCH)
summary (lineoCH)

## check mean distance to nearest trap etc.
summary(LStraps)

## LStraps has several site covariates; terse descriptions are in
##  an extra attribute that may be displayed thus
attr(LStraps, "habitat.variables")

## For density modelling we need covariate values at each point in the
## habitat mask. This requires both on-grid interpolation and
## extrapolation beyond the grids. One (crude) possibility is to
## extrapolate a mask covariate from a covariate of the nearest trap:

LSmask <- make.mask(LStraps, buffer = 30, type = "trapbuffer")
temp <- nearesttrap(LSmask, LStraps)
habclass <- covariates(LStraps)$class[temp]
habclass <- factor (habclass, levels = c(1,2))
covariates(LSmask) <- data.frame(habclass)

## plot mask with colour-coded covariate
par(fg = "white")  ## white pixel borders
plot (LSmask, covariate = "habclass", dots = FALSE, axes = FALSE,
    col = c("yellow", "green"), border = 0)
plot(LStraps, add = TRUE, detpar = list(pch = 16))
par(fg = "black")    ## default
}
