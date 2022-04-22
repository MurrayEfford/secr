\name{extractMoves}
\alias{extractMoves}

\title{Simulated Movements}

\description{

Extract movements from a previously simulated multi-session population.

}

\usage{

extractMoves(pop, plotn = NULL, add = FALSE, collapse = TRUE, maxradius = Inf, ...)

}

\arguments{
  \item{pop}{popn object from \code{\link{sim.popn}}}
  \item{plotn}{integer maximum number of instances to plot at each session}
  \item{add}{logical for whether to add to existing plot}
  \item{collapse}{logical; if TRUE plots for sessions 2, 3,... are added to the first}
  \item{maxradius}{numeric radius for selecting subset of initial locations}
  \item{\dots}{arguments passed to \code{\link{arrows}}}
}

\details{

This function is mostly used to check the movement simulations.

Moves are constrained by the edge (argument `edgemethod' of \code{\link{sim.popn}}). `maxradius' may be set to restrict the extraction to the subset of animals initially near the centroid of the arena in each session.

Plotting uses the graphics function \code{\link{arrows}} that has some quirks, such as difficult-to-suppress warnings for zero-length moves. Set \code{code = 0} to suppress arrowheads; \code{length = 0.1} to shorten to 0.1 inches, etc.

}

\value{

List of data frames, one for each session but the last (columns `x1',`y1',`x2',`y2',`d').

}

\seealso{

\link{turnover}, \code{\link{sim.popn}}

}

\examples{
set.seed(12345)
pop3 <- sim.popn(D = 2, core = make.grid(), buffer = 200, nsessions = 3, 
    details = list(lambda = 1.0, movemodel = 'BVE', move.a = 50, 
    edgemethod = 'stop'))
m <- extractMoves(pop3, plotn = 10, length = 0.1)
mean(unlist(sapply(m, '[', 'd')))    # less than nominal 2 x move.a

# For distances closer to nominal for BVE (2 x move.a = 100), 
# increase size of arena (e.g., buffer = 500) and consider only 
# central animals (e.g., maxradius = 300).

}
