\name{hornedlizard}
\docType{data}

\alias{hornedlizard}
\alias{hornedlizardCH}

\title{ Flat-tailed Horned Lizard Dataset }

\description{

Data from multiple searches for flat-tailed horned lizards
(\emph{Phrynosoma mcalli}) on a plot in Arizona, USA.

}

\usage{hornedlizardCH}

\details{

The flat-tailed horned lizard (\emph{Phrynosoma mcalli}) is a desert
lizard found in parts of southwestern Arizona, southeastern California
and northern Mexico. There is considerable concern about its
conservation status. The species is cryptically coloured and has the
habit of burying under the sand when approached, making it difficult or
impossible to obtain a complete count (Grant and Doherty 2007).
  
K. V. Young conducted a capture--recapture survey of flat-tailed horned
lizards 25 km south of Yuma, Arizona, in the Sonoran Desert. The habitat
was loose sand dominated by creosote bush and occasional bur-sage and
Galletta grass. A 9-ha plot was surveyed 14 times over 17 days (14 June
to 1 July 2005). On each occasion the entire 300 m x 300 m plot was
searched for lizards. Locations within the plot were recorded by
handheld GPS. Lizards were captured by hand and marked individually on
their underside with a permanent marker. Marks are lost when the lizard
sheds, but this happens infrequently and probably caused few or no
identification errors during the 2.5-week study.

A total of 68 individuals were captured 134 times. Exactly half of the
individuals were recaptured at least once.

Royle and Young (2008) analysed the present dataset to demonstrate a
method for density estimation using data augmentation and MCMC
simulation. They noted that the plot size was much larger than has been
suggested as being practical in operational monitoring efforts for this
species, that the plot was chosen specifically because a high density of
individuals was present, and that high densities typically correspond to
less movement in this species. The state space in their analysis was a
square comprising the searched area and a 100-m buffer (J. A. Royle
pers. comm.).

The detector type for these data is `polygonX' and there is a single
detector (the square plot). The data comprise a capture history matrix
(the body of \code{hornedlizardCH}) and the x-y coordinates of each
positive detection (stored as an attribute that may be displayed with
the `xy' function); the `traps' attribute of \code{hornedlizardCH}
contains the vertices of the plot. See
\href{../doc/secr-datainput.pdf}{secr-datainput.pdf} for guidance on
data input.

Non-zero entries in a polygonX capture-history matrix indicate the
number of the polygon containing the detection. In this case there was
just one polygon, so entries are 0 or 1. No animal can appear more than
once per occasion with the polygonX detector type, so there is no need
to specify `binomN = 1' in secr.fit.

\tabular{ll}{
  Object \tab Description \cr
  hornedlizardCH \tab single-session capthist object \cr
}

}

\source{

Royle and Young (2008) and J. A. Royle (pers. comm.), with additional
information from K. V. Young (pers. comm.).

}

\references{

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture analysis of data from area
searches. \emph{Ecology} \bold{92}, 2202--2207.

Grant, T. J. and Doherty, P. F. (2007) Monitoring of the flat-tailed
horned lizard with methods incorporating detection
probability. \emph{Journal of Wildlife Management} \bold{71}, 1050--1056

Marques, T. A., Thomas, L. and Royle, J. A. (2011) A hierarchical model
for spatial capture--recapture data: Comment. \emph{Ecology} \bold{92},
526--528.

Royle, J. A. and Young, K. V. (2008) A hierarchical model for spatial
capture--recapture data. \emph{Ecology} \bold{89}, 2281--2289.

}

\seealso{
\code{\link{capthist}}, \code{\link{detector}}, \code{\link{reduce.capthist}}
}

\examples{

plot(hornedlizardCH, tracks = TRUE, varycol = FALSE,
    lab1 = TRUE, laboff = 6, border = 10, title =
    "Flat-tailed Horned Lizards (Royle & Young 2008)")

table(table(animalID(hornedlizardCH)))
traps(hornedlizardCH)

## show first few x-y coordinates
head(xy(hornedlizardCH))

\dontrun{

## Compare default (Poisson) and binomial models for number
## caught
FTHL.fit <- secr.fit(hornedlizardCH)
FTHLbn.fit <- secr.fit(hornedlizardCH, details =
    list(distribution = "binomial"))
collate(FTHL.fit, FTHLbn.fit)[,,,"D"]

## Collapse occasions (does not run faster)
hornedlizardCH.14 <- reduce(hornedlizardCH, newoccasions =
    list(1:14), outputdetector = "polygon")
FTHL14.fit <- secr.fit(hornedlizardCH.14, binomN = 14)

}

}

\keyword{datasets}
