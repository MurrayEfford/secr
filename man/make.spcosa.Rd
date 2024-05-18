\name{make.spcosa}
\alias{make.spcosa}

\title{Construct Spatial Coverage Design}

\description{

A spatial coverage design places one cluster of detectors in each compact subregion of a region of interest. Equal subregions are determined by \emph{k}-means clustering of pixels (Walvoort et al. 2010).

}

\usage{

make.spcosa(n, cluster, region, rotate = 0, randomize = FALSE, maxtries = 100, keep.mask = FALSE, ...)

}

\arguments{
  \item{n}{integer number of subregions}
  \item{cluster}{traps object defining a cluster of detectors}
  \item{region}{boundary of region of interest}
  \item{rotate}{numeric angular rotation of each cluster (negative for random)}
  \item{randomize}{logical; if TRUE then cluster is located at random within subregion}
  \item{maxtries}{integer maximum attempts to find random location}
  \item{keep.mask}{logical; if TRUE then a habitat mask with covariate `stratum'}
  \item{\dots}{other arguments passed to \code{\link{make.mask}}}
}

\details{

The region may be specified in any form acceptable as the `poly' argument of \code{\link{make.mask}} (see also \link{boundarytoSF}).

The \dots argument determines the coarseness of the discretization used to define the subregions, via the `nx' or `spacing' arguments of \code{\link{make.mask}}. 

By default  (randomize = FALSE) clusters are centred at subregion centroids. Otherwise (randomize = TRUE) clusters are centred in a randomly selected cell of each subregion, subject to the constraint that all detectors fall within the subregion. An error results if no cluster meeting the constraint is found in `maxtries' attempts.

Slightly different partitions of `poly' are generated depending on the value of the random seed, so for consistency this should be first set with \code{set.seed}.

}

\value{

A traps object with n x nrow(cluster) detectors.

}

\references{

Walvoort, D., Brus, D., and de Gruijter, J. (2010) An R package for spatial
coverage sampling and random sampling from compact geographical strata by k-means. 
\emph{Computers & Geosciences} 36:1261--1267.

}

\seealso{

\code{\link{make.mask}},  
\code{\link{boundarytoSF}},  
\link{traps}

}

\examples{

# preliminaries
polygonfile <- system.file("extdata/possumarea.txt", package = "secr")
poly <- read.table(polygonfile, header = TRUE)
subgrid <- make.grid(3,3, spacing = 80)
set.seed(123)

# nx and keep.mask refer to the discretized region of interest
tr <- make.spcosa(n = 5, subgrid, poly, nx = 32, randomize = TRUE, 
          keep.mask = TRUE)

plot(attr(tr,'mask'), dots = FALSE, cov = 'stratum', legend = FALSE)
plot(tr, add = TRUE)
lines(poly)


}
