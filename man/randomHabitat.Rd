\name{randomHabitat}
\alias{randomHabitat}

\title{Random Landscape}

\description{

  The Modified Random Cluster algorithm of Saura and Martinez-Millan
(2000) is used to generate a mask object representing patches of
contiguous `habitat' cells (pixels) within a `non-habitat'
matrix (`non-habitat' cells are optionally dropped). Spatial autocorrelation
(fragmentation) of habitat patches is controlled via the parameter
`p'. `A' is the expected proportion of `habitat' cells.

}

\usage{

randomHabitat(mask, p = 0.5, A = 0.5, directions = 4, minpatch = 1,
drop = TRUE, covname = "habitat", plt = FALSE, seed = NULL)

}

\arguments{
  \item{mask}{secr mask object to use as template}
  \item{p}{parameter to control fragmentation}
  \item{A}{parameter for expected proportion of habitat}
  \item{directions}{integer code for adjacency (rook's move 4 or queen's move 8)}
  \item{minpatch}{integer minimum size of patch}
  \item{drop}{logical for whether to drop non-habitat cells}
  \item{covname}{character name of covariate when \code{drop = FALSE}}
  \item{plt}{logical for whether intermediate stages should be plotted}
  \item{seed}{either NULL or an integer that will be used in a call to \code{set.seed}}
}

\details{

  Habitat is simulated within the region defined by the cells of
  \link{mask}. The region may be non-rectangular.
  
The algorithm comprises stages A-D:

A. Randomly select proportion \code{p} of cells from the input mask

B. Cluster selected cells with any immediate neighbours as defined by
\code{directions}

C. Assign clusters to `non-habitat' (probability 1--A) and `habitat'
(probability A)

D. Cells not in any cluster from (B) receive the habitat class of the
majority of the <=8 adjacent cells assigned in (C), if there are any;
otherwise they are assigned at random (with probabilities 1--A, A).

Fragmentation declines, and cluster size increases, as p increases up to the
`percolation threshold' which is about 0.59 in the default case (Saura
and Martinez-Millan 2000 p.664).

If \code{minpatch > 1} then habitat patches of less than \code{minpatch}
cells are converted to non-habitat, and vice versa. This is likely to
cause the proportion of habitat to deviate from \code{A}.

If \code{drop = FALSE} a binary-valued (0/1) covariate with the
requested name is included in the output mask, which has the same extent
as the input. Otherwise, non-habitat cells are dropped and no covariate
is added.

}

\value{

  An object of class `mask'. By default (\code{drop = TRUE}) this
  has fewer rows (points) than the input mask.
  
  The attribute ``type'' is a character string formed from \code{paste('MRC p=',p, ' A=',A, sep='')}.
  
  The \link{RNG} seed is stored as attribute `seed' (see \link{secrRNG}).

}

\references{

Hijmans, R. J. and van Etten, J. (2011) raster: Geographic analysis and
modeling with raster data. R package version 1.9-33.
\url{https://CRAN.R-project.org/package=raster}.

Saura, S. and Martinez-Millan, J. (2000) Landscape patterns simulation
with a modified random clusters method. \emph{Landscape Ecology},
\bold{15}, 661--678.

}

\note{

  Single-linkage clustering and adjacency operations use functions
  `clump' and `adjacency' of the package \pkg{raster}; `clump' also
  requires package \pkg{igraph0} (\pkg{raster} still uses this
  deprecated version). Optional plotting
  of intermediate stages (\code{plt = TRUE}) uses the plot method for
  rasterLayers in \pkg{raster}.

A non-rectangular input mask is padded out to a rectangular rasterLayer
for operations in \pkg{raster}; cells added as padding are ultimately
dropped.

The procedure of Saura and Martinez-Millan (2000) has been followed as
far as possible, but this implementation may not match theirs in every
detail.

This implementation allows only two habitat classes. The parameter A is
  the \emph{expected} value of the habitat proportion; the
  \emph{realised} habitat proportion may differ quite strongly from A,
  especially for large p (e.g., p > 0.5).

  Anisotropy is not implemented; it would require skewed adjacency
  filters (i.e. other than rook- or queen-move filters) that are not
  available in \pkg{raster}.

Gaussian random fields provide an alternative method for simulating 
random habitats (e.g., LGCP option in \code{\link{sim.popn}}).

}

\seealso{
  
  \code{\link{mask}}, \code{\link{make.mask}}, \code{\link{sim.popn}}
  
}

\examples{

\dontrun{

tempmask <- make.mask(nx = 100, ny = 100, spacing = 20)
mrcmask <- randomHabitat(tempmask, p = 0.4, A = 0.4)
plot(mrcmask, dots = FALSE, col = "green")
pop <- sim.popn(10, mrcmask, model2D = "IHP")
plot(pop, add = TRUE)

## plot intermediate steps A, C, D
opar <- par(mfrow = c(1,3))
mrcmask <- randomHabitat(tempmask, p = 0.4, A = 0.4, plt = TRUE)
par(opar) 

## keep non-habitat cells
mrcmask <- randomHabitat(tempmask, p = 0.4, A = 0.4, drop = FALSE)
plot(mrcmask, covariate = "habitat", dots = FALSE,
    col = c("grey","green"), breaks = 2)

## effect of purging small patches
opar <- par(mfrow=c(1,2))
mrcmask <- randomHabitat(tempmask, p = 0.4, A = 0.4, minpatch = 1)
plot(mrcmask, dots = FALSE, col  ="green")
mrcmask <- randomHabitat(tempmask, p = 0.4, A = 0.4, minpatch = 5)
plot(mrcmask, dots = FALSE, col  ="green")
par(opar)
}

}

\keyword{ datagen }
