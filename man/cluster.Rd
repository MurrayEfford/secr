\name{cluster}
\alias{cluster}
\alias{clusterID}
\alias{clustertrap}
\alias{clusterID<-}
\alias{clustertrap<-}

\title{ Detector Clustering }

\description{
  
  Clusters are uniform groups of detectors. Use these functions to
  extract or replace cluster information of a \code{traps} object, or
  extract cluster information for each detection in a \code{capthist}
  object.
  
}

\usage{
clusterID(object)
clusterID(object) <- value
clustertrap(object)
clustertrap(object) <- value
}
\arguments{
  
  \item{object}{ \code{traps} or \code{capthist} object }
  
  \item{value}{ factor (\code{clusterID}) or integer-valued vector
  (\code{clustertrap}) }

}
\details{

  Easy access to attributes used to define compound designs, those in
  which a detector array comprises several similar subunits
  (`clusters'). `clusterID' identifies the detectors belonging to each
  cluster, and `clustertrap' is a numeric index used to relate matching
  detectors in different clusters.
  
  For replacement (`traps' only), the number of rows of \code{value}
  must match exactly the number of detectors in \code{object}.

  `clusterID' and `clustertrap' are assigned automatically by
  \code{\link{trap.builder}}.

}

\value{

  Factor (\code{clusterID}) or integer-valued vector
  (\code{clustertrap}).

  \code{clusterID(object)} may be NULL.
  
}

\seealso{

\code{\link{traps}}, \code{\link{trap.builder}}, \code{\link{mash}},
\code{\link{derivedCluster}}, \code{\link{cluster.counts}},
\code{\link{cluster.centres}}

}

\examples{

## 25 4-detector clusters
mini <- make.grid(nx = 2, ny = 2)
tempgrid <- trap.builder (cluster = mini , method = "all",
    frame = expand.grid(x = seq(100, 500, 100), y = seq(100,
    500, 100)))
clusterID(tempgrid)
clustertrap(tempgrid)

tempCH <- sim.capthist(tempgrid)
table(clusterID(tempCH)) ## detections per cluster
cluster.counts(tempCH)   ## distinct individuals

}

\keyword{ manip }

