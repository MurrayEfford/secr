\name{make.spcosa}
\alias{make.spcosa}

\title{Construct Spatial Coverage Design}

\description{

A spatial coverage design places one cluster of detectors in each compact subregion of a region of interest. Equal subregions are determined by \emph{k}-means clustering of pixels (Walvoort et al. 2010).

}

\usage{

make.spcosa(n, cluster, poly, rotate = 0, randomize = FALSE, keep.mask = FALSE, ...)

}

\arguments{
  \item{n}{integer number of subregions}
  \item{cluster}{traps object defining a cluster of detectors}
  \item{poly}{boundary of region of interest}
  \item{rotate}{numeric angular rotation of each cluster (negative for random)}
  \item{randomize}{logical; if TRUE then cluster is located at random within subregion (Not yet available)}
  \item{keep.mask}{logical; if TRUE then a habitat mask with covariate `stratum'}
  \item{\dots}{other arguments passed to \code{\link{make.mask}}}
}

\details{

The \dots argument determines the coarseness of the discretization used to define the subregions, via the `nx' or `spacing' arguments of \code{\link{make.mask}}. 

By default, clusters are centred at subregion centroids. 

Otherwise, clusters are placed at uniform random centres within each subregion, subject to the constraint that all detectors fall within the subregion.

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

\link{traps}

}

\examples{

polygonfile <- system.file("extdata/possumarea.txt", package = "secr")
poly <- read.table(polygonfile, header = TRUE)
subgrid <- make.grid(3,3, spacing = 80)
tr <- make.spcosa(5, subgrid, poly, nx = 32, keep.mask = TRUE)
plot(attr(tr,'mask'), dots = FALSE, cov = 'stratum', legend = FALSE)
plot(tr, add = TRUE)


}
