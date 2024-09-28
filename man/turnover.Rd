\name{turnover}
\alias{turnover}

\title{ Specifying a Dynamic Population }

\description{

\code{sim.popn} can simulate a multi-session population with known between-session survival, recruitment and movement probabilities. The parameter settings to achieve this are passed to \code{sim.popn} in its `details' argument. Components of `details' that are relevant to turnover are described below; see \code{\link{sim.popn}} for others.

Multi-session populations are generated in \code{sim.popn} whenever its argument `nsessions' is greater than 1. If \code{details$lambda} remains NULL (the default) then the population for each successive session is generated de novo from the given density model (model2D, D etc.). If a value is specified for \code{details$lambda} then only the first population is generated de novo; remaining populations are generated iteratively with probabilistic mortality, recruitment and movement as described here.

}

\section{Turnover components of \code{sim.popn} details argument}{

\tabular{lll}{
Component \tab Description \tab Default \cr
\code{phi} \tab per capita survival rate \eqn{\phi} \tab 0.7 \cr
\code{survmodel} \tab probability model for number of survivors \tab ``binomial" \cr
\code{lambda} \tab finite rate of increase \eqn{\lambda = N_{t+1} / N_t} \tab none \cr
\code{recrmodel} \tab probability model for number of recruits \tab ``poisson" \cr
\code{superN} \tab optional superpopulation size for `multinomial' recruitment model\tab see below \cr
\code{Nrecruits} \tab number of recruits to add at t+1 for `specified' recruitment model \tab 0 \cr
\code{movemodel} \tab ``static", ``uncorrelated", ``normal", ``exponential", ``t2D" or a user function \tab ``static" \cr
\code{move.a} \tab first parameter of movement kernel (replacing sigma.m) \tab 0 \cr
\code{move.b} \tab second parameter of movement kernel \tab 1 \cr

\code{edgemethod} \tab treatment of animals that cross the boundary \tab ``wrap" \cr
\code{sigma.m} \tab deprecated in 3.2.1; use move.a \tab 0 \cr
\code{wrap} \tab deprecated in 3.1.6; use edgemethod \tab TRUE i.e. edgemethod = ``wrap" \cr
}

}

\section{Survival}{

Survival is usually thought of as a Bernoulli process (outcome 0 or 1 for each individual) so the number of survivors \eqn{S} is a binomial variable (survmodel = "binomial"). Another approach is to fix the proportion surviving, but this can be done exactly only when \eqn{\phi N} is an integer. A (slightly ad hoc) solution is to randomly choose between the two nearest integers with probability designed in the long term (over many sessions) to give the required \eqn{\phi} (survmodel = "discrete").

}

\section{Population growth and recruitment}{

Per capita recruitment (\eqn{f}) is the difference between lambda and phi (\eqn{f = \lambda - \phi}), which must be non-negative (phi > lambda causes an error). The number of recruits B is a random variable whose probability distribution is controlled by details$recrmodel:
\tabular{ll}{
Value \tab Probability model \cr
"constantN" \tab Exact replacement of animals that die (B = \eqn{N_t - S})\cr
"binomial" \tab Binomial number of recruits (B ~ bin(\eqn{N_t, f}) \cr
"poisson" \tab Poisson number of recruits (B ~ pois(\eqn{f N_t})) \cr
"discrete" \tab Minimum-variance number of recruits (see Survival) \cr
"multinomial" \tab The POPAN model, conditioning on superpopulation size (e.g., Schwarz and Arnason 1996)) \cr
"specified" \tab Add the number of recruits specified in Nrecruits (may be vector) \cr
}

In the case of binomial recruitment there is a maximum of one recruit per existing individual, so lambda <= (phi+1). Multinomial recruitment requires a value for the superpopulation size. This may be provided as the details component "superN". If not specified directly, a value is inferred by projecting a trial initial (simulated) population using the specified phi and lambda.

Specifying the integer number of recruits in each year (recrmodel `specified') overrides the value of lambda, but a non-null value should be given for lambda.

}

\section{Movement}{

Individuals may shift their home range centre between sessions. Movement probability is governed by a circular kernel specified by `movemodel' and the parameter values `move.a' and `move.b' (optional). By default there is no movement between sessions (movemodel = "static"). Other options are

\tabular{lll}{
``IND" \tab ``uncorrelated" \tab individuals are randomly assigned a new, independent location within the buffered area\cr
``BVN" \tab ``normal"\tab bivariate normal (Gaussian) kernel with parameter move.a (previously called sigma.m)\cr
``BVE" \tab ``exponential" \tab negative exponential (Laplace) kernel with parameter move.a \cr
``BVT" \tab ``t2D"\tab circular 2-dimensional t-distribution with scale parameter move.a and shape parameter move.b = df/2 (2Dt of Clark et al. 1999)\cr
``RDE" \tab \tab exponential distribution of radial distance (Ergon & Gardner, 2014) \cr
``RDG" \tab \tab gamma distribution of radial distance (Ergon & Gardner, 2014) \cr
``RDL" \tab \tab log-normal distribution of radial distance (Ergon & Gardner, 2014) \cr
\tab \tab  (parameterized with move.a = exp(mu), move.b = 1/CV^2 = 1 / (exp(SD^2) - 1) \cr

}

The package \pkg{openCR} >=1.4.0 provides functions for constructing and plotting these kernels and summarising their properties (\code{make.kernel}; \code{plot} and \code{summary} methods for kernel objects). The \pkg{secr} function \code{\link{extractMoves}} is useful for checking simulations of movement.

Models IND, BVN, BVE, and RDE may also be zero-inflated (suffix ``zi"). The parameter `move.a' (INDzi) or `move.b' (BVNzi, BVEzi, RDEzi) is the zero-inflation probability. See Examples.


In \pkg{secr} <3.2.1 sigma.m was also used to indicate two special cases; these continue to work but may be discontinued in the future:

sigma.m = 0 corresponds to movemodel = `static'

sigma.m < 0 corresponds to movemodel = `uncorrelated'

In \pkg{secr} >= 4.4.0, the `movemodel' component may also be a user-provided function with these characteristics: two or three arguments, the first being the number of centres to be moved (e.g., n) and the others parameters of the dispersal distribution (e.g., a,b); the function should return a matrix of n rows and 2 columns, the displacements in the x- and y-directions. The output is a set of random points from the bivariate dispersal kernel. The function will be called with the current number of centres and parameter values move.a and move.b as needed.

If movement takes an animal across the boundary of the arena (buffered area) in \code{sim.popn} the component "edgemethod" comes into play. By default, locations are toroidally wrapped i.e. the animal re-joins the population on the opposing edge. Other options are ``clip'' (discard), ``clipandreplace'' (assign new identity at original location), ``stop'' (stop just inside the boundary), ``reflect'' (bounce off edges to the limit of the dispersal), ``normalize'' = ``truncate'' (truncate kernel and scale probability to 1.0) and ``none" (allow centres outside the buffered area). The ``normalize'' option (new in \pkg{secr} 4.3.3) can take longer as it repeatedly relocates each individual until its destination lies within the bounding box, up to a maximum of 500 attempts.

}

\seealso{
  \code{\link{sim.popn}}, \code{\link{extractMoves}}
}

\references{

Clark, J. S, Silman, M., Kern, R., Macklin, E. and HilleRisLambers, J. (1999) Seed dispersal near and far: patterns across temperate and tropical forests. \emph{Ecology} \bold{80}, 1475--1494.

Nathan , R., Klein, E., Robledo-Arnuncio, J. J. and Revilla, E. (2012) 
Dispersal kernels: a review. In: J Clobert et al. \emph{Dispersal Ecology and Evolution}. Oxford University Press. Pp 187--210.

}

\examples{

par (mfrow = c(2,3), mar = c(1,1,1,1))

## birth and death only
grid <- make.grid(nx = 7, ny = 4, detector = 'proximity', spacing = 10)
pop <- sim.popn (Nbuffer = 100, core = grid, nsessions = 6,    
    details = list(lambda = 0.8, phi = 0.6))
sapply(pop, nrow)  ## how many individuals?
plot(pop)

## movement only
pop2 <- sim.popn (Nbuffer = 100, core = grid, nsessions = 6,    
    details = list(lambda = 1, phi = 1, movemodel = 'normal', 
    move.a = 10, edgemethod = "wrap"))
pop3 <- sim.popn (Nbuffer = 100, core = grid, nsessions = 6,    
    details = list(lambda = 1, phi = 1, movemodel = 'normal', 
    move.a = 10, edgemethod = "clip"))
pop4 <- sim.popn (Nbuffer = 100, core = grid, nsessions = 10,    
    details = list(lambda = 1, phi = 1, movemodel = 'normal', 
    move.a = 10, edgemethod = "stop"))
sapply(pop2, nrow)  ## how many individuals?
plot(pop2)

## show effect of edgemethod --
## first session blue, last session red
cols <- c('blue',rep('white',4),'red')
par (mfrow=c(1,2))
plot(pop2, collapse = TRUE, seqcol = cols)
plot(pop3, collapse = TRUE, seqcol = cols)

## zero-inflated movement
## move.b is zero-inflation probability
pop5 <- sim.popn (Nbuffer = 1000, core = grid, nsessions = 6,    
     details = list(lambda = 1, phi = 1, movemodel = 'RDEzi', 
         move.a = 50, move.b = 0.5, edgemethod = "none"))
mean(do.call(rbind,extractMoves(pop5))$d)   # approx 50 * 0.5

}
