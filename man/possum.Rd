\name{possum}
\docType{data}

\alias{possum}
\alias{possumCH}
\alias{possumarea}
\alias{possumremovalarea}
\alias{possummask}
\alias{possum.model.0}
\alias{possum.model.Ds}

\title{ Brushtail Possum Trapping Dataset }

\description{

Data from a trapping study of brushtail possums at Waitarere, North
Island, New Zealand.

}

\usage{

possumCH
possumarea
possumremovalarea
possummask
possum.model.0
possum.model.Ds

}

\details{

Brushtail possums (\emph{Trichosurus vulpecula}) are an unwanted
invasive species in New Zealand.  Although most abundant in forests,
where they occasionally exceed densities of 15 / ha, possums live
wherever there are palatable food plants and shelter.

Efford et al. (2005) reported a live-trapping study of possums in
\emph{Pinus radiata} plantation on coastal sand dunes.  The 300-ha site
at Waitarere in the North Island of New Zealand was a peninsula, bounded
on one side by the sea and on two other sides by the Manawatu river. 
Cage traps were set in groups of 36 at 20-m spacing around the perimeter
of five squares, each 180 m on a side.  The squares (`hollow grids') were
centred at random points within the 300-ha area.  Animals were tagged and
released daily for 5 days in April 2002.  Subsequently, leg-hold trapping
was conducted on a trapping web centred on each square (data not
reported here), and strenuous efforts were made to remove all possums by
cyanide poisoning and further leghold trapping across the entire area. 
This yielded a density estimate of 2.26 possums / ha.

Traps could catch at most one animal per day.  The live-trapped animals
comprised 46 adult females, 33 adult males, 10 immature females and 11
immature males; sex and/or age were not recorded for 4 individuals (M.
Coleman unpubl. data). These counts do not sum to the number of capture
histories - see Note.  One female possum was twice captured at two sites
on one day, having entered a second trap after being released; one
record in each pair was selected arbitrarily and discarded.

The data are provided as a single-session \code{capthist} object
`possumCH'.  `possummask' is a matching mask object - see Examples.
Two fitted models are provided for illustration.

The dataframe \code{possumarea} contains boundary coordinates of a
habitat polygon that is used to clip \code{possummask} at the shore
(from secr 1.5). \code{possumarea} comprises a single polygon
representing the extent of terrestrial vegetation to the west, north and
east, and an arbitrary straight southern boundary. The boundary is also
included as a shapefile and as a text file (`possumarea.shp' etc. and
`possumarea.txt' in the package `extdata' folder). See Examples in
\code{\link{make.mask}}.

The dataframe \code{possumremovalarea} contains boundary coordinates of
another polygon, the nominal removal area of Efford et al. (2005 Fig. 1)
(from secr 2.3).

\tabular{ll}{
  Object \tab Description \cr
  possumCH \tab capthist object \cr
  possummask \tab mask object \cr
  possumarea \tab habitat perimeter \cr
  possumremovalarea \tab nominal boundary of removal region \cr
  possum.model.0 \tab fitted secr model -- null \cr
  possum.model.Ds \tab fitted secr model -- distance to shore \cr
}

}

\source{

Landcare Research, New Zealand.

}

\note{

  A significant problem with the data used by Efford et al. (2005) was
noticed recently. Five capture histories in possumCH are for animals
that had lost a previous tag. A further three histories may also have
been animals that were tagged previously or mis-recorded. Analyses that
treat each previously tagged animal as a new individual are in error
(this includes the published analyses, the pre-fitted models described
here, and those in the vignette secr-densitysurfaces.pdf). All eight
questionable histories are now indicated in possumCH with the logical
covariate `prev.tagged'.

Methods have not yet been developed to adjust for tag loss in SECR
models. 
}
  
\references{

Borchers, D.L. and Efford, M.G. (2008) Spatially explicit maximum
likelihood methods for capture-recapture studies. \emph{Biometrics} \bold{64},
377--385.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
\emph{Animal Biodiversity and Conservation} \bold{27}, 217--228. 

Efford, M. G., Warburton, B., Coleman, M. C. and Barker, R. J. (2005) A
field test of two methods for density estimation. \emph{Wildlife Society
Bulletin} \bold{33}, 731--738.

}

\seealso{

\code{\link{capthist}}

}

\examples{

plot(possummask)
plot(possumCH, tracks = TRUE, add = TRUE)
plot(traps(possumCH), add = TRUE)
lines(possumarea)
summary(possumCH)

## compare & average pre-fitted models
AIC(possum.model.0, possum.model.Ds)
modelAverage(possum.model.0, possum.model.Ds)

\dontrun{

## Roughly estimate tag-loss error by dropping dubious histories
## i.e. restrict to "not previously tagged"
NPT <- !covariates(possumCH)$prev.tagged
possum.model.0.NPT <- secr.fit(subset(possumCH,NPT), mask =
    possummask, trace = FALSE)
predict(possum.model.0)[1,2]/ predict(possum.model.0.NPT)[1,2]
## ...about 9\%
}

}

