\name{trap.builder}
\alias{trap.builder}
\alias{mash}
\alias{cluster.counts}
\alias{cluster.centres}

\title{
Complex Detector Layouts
}

\description{
  
Construct detector layouts comprising small arrays (clusters) replicated
across space, possibly at a probability sample of points.

}

\usage{

trap.builder (n = 10, cluster, region = NULL, frame = NULL, 
    method = c("SRS", "GRTS", "all", "rank"), 
    edgemethod = c("clip", "allowoverlap", "allinside", "anyinside", "centreinside"), 
    samplefactor = 2, ranks = NULL, rotation = NULL, detector,
    exclude = NULL, exclmethod = c("clip", "alloutside", "anyoutside", "centreoutside"), 
    plt = FALSE, add = FALSE, \dots)

mash (object, origin = c(0,0), clustergroup = NULL, \dots)

cluster.counts (object)

cluster.centres (object)

}

\arguments{
  \item{n}{ integer number of clusters (ignored if method = "all") }
  \item{cluster}{ traps object }
  \item{region}{ bounding polygon(s)}
  \item{frame}{ data frame of points used as a finite sampling frame }
  \item{method}{ character string (see Details) }
  \item{edgemethod}{character string (see Details) }
  \item{samplefactor}{ oversampling to allow for rejection of edge clusters
    (multiple of n) }
  \item{ranks}{ vector of relative importance (see Details)}
  \item{rotation}{ angular rotation of each cluster about centre (degrees) }
  \item{detector}{ character detector type (see \code{\link{detector}})}
  \item{exclude}{ polygon(s) from which detectors are to be excluded}
  \item{exclmethod}{character string (see Details)}
  \item{plt}{ logical: should array be plotted? }
  \item{add}{ logical: add to existing plot }
  \item{object}{ single-session multi-cluster capthist object, or traps
    object for \code{cluster.centres} }
  \item{origin}{ new coordinate origin for detector array }
  \item{clustergroup}{ list of vectors subscripting the clusters to be
    mashed }
  \item{\dots}{other arguments passed by trap.builder to spsurvey::grts (e.g., mindis) 
  and by mash to make.capthist (e.g., sortrows)}
}

\details{

  The detector array in \code{cluster} is replicated \code{n}
  times and translated to centres sampled from the area sampling frame
  in \code{region} or the finite sampling frame in \code{frame}. Each
  cluster may be rotated about its centre either by a fixed number of
  degrees (\code{rotation} positive), or by a random angle (\code{rotation}
  negative).
  
  If the \code{cluster} argument is not provided then single detectors of
  the given type are placed according to the design.

  The sampling frame is finite (the points in \code{frame}) whenever
  \code{frame} is not NULL. If \code{region} and \code{frame} are both
  specified, sampling uses the finite frame but sites may be clipped
  using the polygon.

  \code{region} and \code{exclude} may be a two-column matrix or
  dataframe of x-y coordinates for the boundary, or one of the other polygon 
  sources listed in \code{\link{boundarytoSF}} (these allow multiple polygons). 
  
  \code{method} may be "SRS", "GRTS", "all" or "rank". "SRS" takes a simple
  random sample (without replacement in the case of a finite sampling
  frame). "GRTS" takes a spatially representative sample using the
  `generalized random tessellation stratified' (GRTS) method of Stevens
  and Olsen (2004). "all" replicates \code{cluster} across all points in
  the finite sampling frame. "rank" selects \code{n} sites from
  \code{frame} on the basis of their ranking on the vector `ranks',
  which should have length equal to the number of rows in
  \code{frame}; ties are resolved by drawing a site at random.
  
  Options for \code{edgemethod} are --
  \tabular{llll}{
  \code{edgemethod} \tab Description\cr
  \code{"clip"} \tab reject any individual detectors outside \code{region} \cr
  \code{"allowoverlap"} \tab no action \cr
  \code{"allinside"} \tab reject whole cluster if any component is outside \code{region} \cr
  \code{"anyinside"} \tab reject whole cluster if no component is inside \code{region}\cr
  \code{"centreinside"} \tab reject whole cluster if centre outside \code{region}, and clip to \code{region}
  }

  Similarly, \code{exclmethod} may be "clip" (reject individual detectors), 
  "alloutside" (reject whole cluster if any component is outside \code{exclude}) etc. 
  Sufficient additional samples (\code{(samplefactor--1) * n}) must be drawn to 
  allow for replacement of any rejected clusters; otherwise, an error is reported 
  (`not enough clusters within polygon').
  
  GRTS samples require function \code{grts} in version >= 5.3.0 of package \pkg{spsurvey} 
  (Dumelle et al. 2022). More sophisticated stratified designs may be specified by using \code{grts}
  directly.

  \code{mash} collapses a multi-cluster capthist object as if all
  detections were made on a single cluster. The new detector coordinates
  in the `traps' attribute are for a single cluster with (min(x),
  min(y)) given by \code{origin}. \code{clustergroup} optionally selects
  one or more groups of clusters to mash; if \code{length(clustergroup)
  > 1} then a multisession capthist object will be generated, one
`session' per clustergroup. By default, all clusters are mashed.

  \code{mash} discards detector-level covariates and occasion-specific
  `usage', with a warning.

  \code{cluster.counts} returns the number of \emph{distinct}
  individuals detected per cluster in a single-session multi-cluster
  capthist object.

  \code{cluster.centres} returns the centroid of the detector locations in each cluster. When clusters have been truncated these differ from the attribute \code{centres} set by \code{\link{make.systematic}}.

}

\value{

  \code{trap.builder} produces an object of class `traps'.

  \code{plt = TRUE} causes a plot to be displayed, including the polygon
  or finite sampling frame as appropriate.

  \code{mash} produces a capthist object with the same number of rows as
  the input but different detector numbering and `traps'. An attribute
  `n.mash' is a vector of the numbers recorded at each cluster; its
  length is the number of clusters. An attribute `centres' is a
  dataframe containing the x-y coordinates of the cluster centres. The
  \code{predict} method for secr objects and the function \code{derived}
  both recognise and adjust for mashing.

  \code{cluster.counts} returns a vector with the number of individuals
  detected at each cluster.

  \code{cluster.centres} returns a dataframe of x- and y-coordinates.

}

\references{
  
  Dumelle, M., Kincaid, T. M., Olsen, A. R., and Weber, M. H. (2021). 
  spsurvey: Spatial Sampling Design and Analysis. R package version 5.2.0.

  Stevens, D. L., Jr., and Olsen, A. R. (2004) Spatially-balanced
  sampling of natural resources. \emph{Journal of the American
  Statistical Association} \bold{99}, 262--278.

}

\note{

  The function \code{\link{make.systematic}} should be used to generate
  systematic random layouts. It calls \code{trap.builder}.

  The sequence number of the cluster to which each detector belongs, and
  its within-cluster sequence number, may be retrieved with the
  functions \code{\link{clusterID}} and \code{\link{clustertrap}}.

}

\seealso{
\code{\link{make.grid}}, \code{\link{traps}},
\code{\link{make.systematic}},
\code{\link{clusterID}},
\code{\link{clustertrap}}
}

\examples{

## solitary detectors placed randomly within a rectangle
tempgrid <- trap.builder (n = 10, method = "SRS",
    region = cbind(x = c(0,1000,1000,0),
    y = c(0,0,1000,1000)), plt = TRUE)

## one detector in each 100-m grid cell -
## a form of stratified simple random sample
## see also Examples in ?make.grid

origins <- expand.grid(x = seq(0, 900, 100),
    y = seq(0, 1100, 100))
XY <- origins + runif(10 * 12 * 2) * 100
temp <- trap.builder (frame = XY, method = "all",
    detector = "multi")
## same as temp <- read.traps(data = XY)
plot(temp, border = 0)  ## default grid is 100 m

\dontrun{

## simulate some data
## regular lattice of mini-arrays
minigrid <- make.grid(nx = 3, ny = 3, spacing = 50,
    detector = "proximity")
tempgrid <- trap.builder (cluster = minigrid , method =
    "all", frame = expand.grid(x = seq(1000, 5000, 2000),
    y = seq(1000, 5000, 2000)), plt = TRUE)
tempcapt <- sim.capthist(tempgrid, popn = list(D = 10))
cluster.counts(tempcapt)
cluster.centres(tempgrid)

## "mash" the CH
summary(mash(tempcapt))

## compare timings (estimates are near identical)
tempmask1 <- make.mask(tempgrid, type = "clusterrect",
    buffer = 200, spacing = 10)
fit1 <- secr.fit(tempcapt, mask = tempmask1, trace = FALSE)        

tempmask2 <- make.mask(minigrid, spacing = 10)
fit2 <- secr.fit(mash(tempcapt), mask = tempmask2, trace = FALSE)  
## density estimate is adjusted automatically
## for the number of mashed clusters (9)

predict(fit1)
predict(fit2)
fit1$proctime
fit2$proctime

## SRS excluding detectors from a polygon

region <- cbind(x = c(0,6000,6000,0,0), y = c(0,0,6000,6000,0))
exclude <- cbind(x = c(3000,7000,7000,3000,3000), y = c(2000,2000,4000,4000,2000))
newgrid <- trap.builder (n = 40, cluster = minigrid,
    method = "SRS", edgemethod = "allinside", region = region,
    exclude = exclude, exclmethod = "alloutside",
    plt = TRUE)

## two-phase design: preliminary sample across region,
## followed by selection of sites for intensive grids

arena <- data.frame(x = c(0,2000,2000,0), y = c(0,0,2500,2500))
t1 <- make.grid(nx = 1, ny = 1)
t4 <- make.grid(nx = 4, ny = 4, spacing = 50)
singletraps <- make.systematic (n = c(8,10), cluster = t1,
    region = arena)
CH <- sim.capthist(singletraps, popn = list(D = 2))
plot(CH, type = "n.per.cluster", title = "Number per cluster")
temp <- trap.builder(10, frame = traps(CH), cluster = t4,
    ranks = cluster.counts(CH), method = "rank",
    edgemethod = "allowoverlap", plt = TRUE, add = TRUE)

## GRTS sample of mini-grids within a rectangle
## GRTS uses package 'spsurvey' >= 5.3.0

minigrid <- make.grid(nx = 3, ny = 3, spacing = 50,
    detector = "proximity")
region <- cbind(x = c(0,6000,6000,0,0), y = c(0,0,6000,6000,0))

if (requireNamespace("spsurvey", versionCheck = list(version = ">=5.3.0"))) {

  tempgrid <- trap.builder (n = 20, cluster = minigrid, region = region, 
      plt = TRUE, method = "GRTS")
  
  # specifying minimum distance between cluster origins
  tempgrid2 <- trap.builder (n = 20, cluster = minigrid, region = region, 
      plt = TRUE, method = "GRTS", mindis = 500, maxtry = 10)
  # use spsurvey::warnprnt() to view warnings (e.g., maxtry inadequate)
  
}

}

}

