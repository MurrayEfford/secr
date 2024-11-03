\name{OVpossum}
\alias{OVpossum}
\alias{OVpossumCH}
\docType{data}

\title{Orongorongo Valley Brushtail Possums}

\description{

A dataset from long-term capture-recapture trapping of brushtail possums
\emph{Trichosurus vulpecula} in New Zealand.

}

\usage{OVpossumCH}

\format{

  A multi-session capthist object of 6 sessions. Sessions are labeled
  49--54, corresponding to February 1996, June 1996, September 1996,
  February 1997, June 1997 and September 1997.
    
}

\details{

Brushtail possums are 2-4 kg largely arboreal marsupials that have
become pests of forests and farmland in New Zealand since their
introduction from Australia in the nineteenth century. Their population
dynamics in mixed native forest have been studied by capture-recapture
in the Orongorongo Valley near Wellington since 1966 (e.g. Crawley 1973,
Efford 1998, Efford and Cowan 2004).

From 1996 to 2006, a grid of 167 traps set on the ground at 30-m spacing
was operated in an area of about 14 ha for 5 consecutive days three
times each year (Efford and Cowan 2004). Each trap could catch only one
animal, with rare exceptions when a young `backrider' entered the trap
with its mother. All animals were tagged and tattooed for individual
identification and released at the site of capture.

A broad shingle riverbed forms a natural boundary on two sides of the
study grid. Much of the grid lies on a gently sloping old alluvial fan
and recent terraces, but to the southeast the valley side rises steeply
and, except where cut by streams, supports beech forest
(\emph{Nothofagus truncata} and \emph{Nothofagus solandri solandri})
rather than the mixed broadleaf forest of the valley floor.

This dataset relates to six five-day trapping sessions in 1996 and 1997,
a time of high and declining density. Possums are long-lived (up to
about 15 years) and as adults restrict their movements to a home range
of 1-10 ha. Breeding is seasonal, resulting in an influx of newly
independent juveniles in the first trapping of each calendar year.

The dataset includes individual covariates not provided by Efford (2012):
`sex' (F or M) and `Ageclass' (1 for first year, 2 for older).

A coarse habitat map is provided for the immediate vicinity of the
trapping grid as the shapefile `OVforest.shp' in the package `extdata'
folder. This distinguishes two forest classes (`beech' and `non-beech'),
and leaves out the shingle riverbed. The distinction between `beech' and
`non-beech' is mapped only to a distance of about 120 m from the
outermost traps. A text file 'leftbank.txt' in the same folder contains the x- and y- coordinates of the adjoining bank of the Orongorongo River. All coordinates relate to the old New Zealand Map Grid (NZMG), since replaced by the New Zealand Transverse Mercator grid (NZTM2000).

The example code shows how to import the shapefile as a \pkg{sp}
SpatialPolygonsDataFrame object and use it to construct a mask for
\code{secr.fit}.

}

\source{

Efford (2012) and unpublished data.

}

\references{

Crawley, M. C. (1973) A live-trapping study of Australian brush-tailed
     possums, \emph{Trichosurus vulpecula} (Kerr), in the Orongorongo
     Valley, Wellington, New Zealand. \emph{Australian Journal of
     Zoology} \bold{21}, 75--90.

Efford, M. G. (1998) Demographic consequences of sex-biased dispersal in
     a population of brushtail possums. \emph{Journal of Animal Ecology}
     \bold{67}, 503--517.

Efford, M. G. (2012) DENSITY 5.0: software for spatially explicit
     capture-recapture. Department of Mathematics and Statistics,
     University of Otago, Dunedin, New
     Zealand. \url{https://www.otago.ac.nz/density/}

Efford, M. G. and Cowan, P. E. (2004) Long-term population trend of
     \emph{Trichosurus vulpecula} in the Orongorongo Valley, New
     Zealand. In: \emph{The Biology of Australian Possums and
     Gliders}. Edited by R. L. Goldingay and S. M. Jackson. Surrey
     Beatty & Sons, Chipping Norton. Pp. 471--483.
     
Ward, G. D. (1978) Habitat use and home range of radio-tagged opossums 
    \emph{Trichosurus vulpecula} (Kerr) in New Zealand lowland forest. 
    In: \emph{The ecology of arboreal folivores}. Edited by G. G. Montgomery. 
    Smithsonian Institute Press. Washington, D.C. Pp. 267--287.

}

\examples{

\dontrun{

library(sf)

summary(OVpossumCH, terse = TRUE)
ovtrap <- traps(OVpossumCH[[1]])

## retrieve and plot the forest map
OVforest <- st_read(system.file("extdata/OVforest.shp", package = "secr"))
## omit forest across the river by selecting only 1,2
OVforest <- OVforest[1:2,]
forestcol <- terrain.colors(6)[c(4,2)]
plot(st_as_sfc(OVforest), col = forestcol)
plot(ovtrap, add = TRUE)

## construct a mask
ovmask <- make.mask(ovtrap, buffer = 120, type = 'trapbuffer',
    poly = OVforest, spacing = 7.5, keep.poly = FALSE)
ovmask <- addCovariates(ovmask, OVforest)

## display mask
par(mar = c(0,0,0,8))
plot(ovmask, covariate = 'forest', dots = FALSE, col = forestcol)
plot(ovtrap, add = TRUE)

## add the left bank of the Orongorongo River
lines(read.table(system.file("extdata/leftbank.txt", package = "secr")))

}

}

