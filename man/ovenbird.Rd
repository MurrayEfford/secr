\name{ovenbird}
\docType{data}

\alias{ovenbird}
\alias{ovenCH}
\alias{ovenCHp}
\alias{ovenbird.model.1}
\alias{ovenbird.model.D}
\alias{ovenmask}

\title{ Ovenbird Mist-netting Dataset }

\description{

Data from a multi-year mist-netting study of ovenbirds (\emph{Seiurus
aurocapilla}) at a site in Maryland, USA.
}

\usage{

ovenCH
ovenCHp
ovenbird.model.1
ovenbird.model.D
ovenmask

}

\details{

From 2005 to 2009 D. K. Dawson and M. G. Efford conducted a
capture--recapture survey of breeding birds in deciduous forest at the
Patuxent Research Refuge near Laurel, Maryland, USA.  The forest was
described by Stamm, Davis & Robbins (1960), and has changed little
since.  Analyses of data from previous mist-netting at the site by Chan
Robbins were described in Efford, Dawson & Robbins (2004) and Borchers &
Efford (2008).

Forty-four mist nets (12 m long, 30-mm mesh) spaced 30 m apart on the
perimeter of a 600-m x 100-m rectangle were operated for approximately 9
hours on each of 9 or 10 non-consecutive days during late May and June
in each year.  Netting was passive (i.e. song playback was not used to
lure birds into the nets).  Birds received individually numbered bands,
and both newly banded and previously banded birds were released at the
net where captured.  Sex was determined in the hand from the presence of
a brood patch (females) or cloacal protuberance (males).  A small amount
of extra netting was done by other researchers after the main session in
some years.

This dataset comprises all records of adult (after-hatch-year) ovenbirds
caught during the main session in each of the five years 2005--2009.  One
ovenbird was killed by a predator in the net in 2009, as indicated by a
negative net number in the dataset.  Sex was determined in the hand from
the presence of a brood patch (females) or cloacal protuberance (males).
 Birds are listed by their band number (4-digit prefix, `.', and 5-digit
number).  


The data are provided as a multi-session \code{capthist} object
`ovenCHp'.  Sex is coded as a categorical individual covariate ("M"
or "F").

Recaptures at the same site within a day are not included in this dataset, 
so \code{ovenCHp} has detector type `proximity'. Previous versions of \pkg{secr}
provided only a trimmed version of these data, retaining only one capture 
per bird per day (\code{ovenCH} with detector type `multi'). That may be 
obtained from \code{ovenCHp} as shown in the examples.

Although several individuals were captured in
more than one year, no use is made of this information in the analyses
presently offered in \pkg{secr}.


An analysis of the data for males in the first four years showed that
they tended to avoid nets after their first capture within a season
(Dawson & Efford 2009).  While the species was present consistently,
the number of detections in any one year was too small to give reliable
estimates of density; pooling of detection parameters across years
helped to improve precision.

Included with the data are a mask and two models fitted to \code{ovenCH} as in
Examples. 

\tabular{ll}{
  Object \tab Description \cr
  ovenCH \tab multi-session capthist object (as multi-catch) \cr
  ovenCHp \tab multi-session capthist object (as binary proximity) \cr
  ovenbird.model.1 \tab fitted secr model -- null \cr
  ovenbird.model.D \tab fitted secr model -- trend in density across years
  \cr
  ovenmask \tab mask object \cr
}

}

\source{

D. K. Dawson (\email{ddawson@usgs.gov}) and M. G. Efford unpublished data.

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture-recapture studies. \emph{Biometrics} \bold{64},
377--385.

Dawson, D. K. and Efford, M. G. (2009) Bird population density estimated from
acoustic signals. \emph{Journal of Applied Ecology} \bold{46}, 1201--1209.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
\emph{Animal Biodiversity and Conservation} \bold{27}, 217--228.

Stamm, D. D., Davis, D. E. and Robbins, C. S. (1960) A method of studying
wild bird populations by mist-netting and banding. \emph{Bird-Banding} \bold{31},
115--130.

}

\seealso{

\code{\link{capthist}}

}

\examples{

## commands used to create ovenCH from the input files
## "netsites0509.txt" and "ovencapt.txt"
## for information only - these files not distributed

# netsites0509 <- read.traps(file = "netsites0509.txt",
#    skip = 1, detector = "proximity")
# temp <- read.table("ovencapt.txt", colClasses=c("character",
#     "character", "numeric", "numeric", "character"))
# ovenCHp <- make.capthist(temp, netsites0509, covnames = "Sex")
# ovenCHp <- reduce(ovenCHp, dropunused = FALSE) # drop repeat detections

par(mfrow = c(1,5), mar = c(1,1,4,1))
plot(ovenCHp, tracks = TRUE, varycol = TRUE)
par(mfrow = c(1,1), mar = c(5,4,4,2) + 0.1)  ## defaults

counts(ovenCHp, "n")

\dontrun{

## trimmed version of data - for consistency with earlier versions

ovenCH <- reduce(ovenCHp, outputdetector = "multi", dropunused = FALSE)

## array constant over years, so build mask only once
ovenmask <- make.mask(traps(ovenCH)[["2005"]], type = "pdot", 
    buffer = 400, spacing = 15, detectpar = list(g0 = 0.03, 
    sigma = 90), nocc = 10)

## fit constant-density model
ovenbird.model.1 <- secr.fit(ovenCH, mask = ovenmask)

## fit temporal trend in density (Session capitalized)
ovenbird.model.D <- secr.fit(ovenCH, mask = ovenmask, 
    model = list(D ~ Session))

## compare pre-fitted models
AIC(ovenbird.model.1, ovenbird.model.D)

}


}

\keyword{datasets}

%\figure{ovenbird.jpg}{options: width=1.5in}
