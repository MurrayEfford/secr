\name{secrRNG}
\alias{random numbers}
\alias{seed}
\alias{secrRNG}

\title{ Random Number Seed }

\description{

The use of random number seeds in \pkg{secr} is explained.

}

\section{Random numbers in \R}{

\R provides several kinds of random number generator (RNG) in the \pkg{base} package (see \link{RNG}). These are used both explicitly, in functions such as \code{\link{runif}} and \code{\link{rnorm}}, and implicitly (\code{\link{sample}}).

A seed suitable for any kind of RNG is held in a vector of 626 integers named \code{.Random.seed}. The vector is not to be modified directly by users. Instead, to start a reproducible stream of random numbers, the user calls \code{\link{set.seed}} with a single non-null integer argument. This has the effect of initialising \code{.Random.seed}. The value of \code{.Random.seed} may nevertheless be stored and restored to reset the RNG state.

\code{\link{set.seed}} with a NULL argument initialises \code{.Random.seed} to an indeterminate (time- and process-dependent) value. The same happens if a random number function is called before \code{.Random.seed} has been set.

}

\section{Handling of RNG seed for simulation in package \pkg{stats}}{

The `official' approach to setting and storing the RNG seed is shown in code and documentation for the generic function \code{\link{simulate}} in the \pkg{stats} package. 

\itemize{

\item The generic has argument `seed' with default NULL. 

\item If `seed' is non-null then \code{set.seed} is called.

\item The returned value has an attribute ``seed'' whose value is either (i) if specified, the integer value of the `seed' argument (with its own attribute ``kind'' from RNGkind), or (ii) the original vector \code{.Random.seed}.

\item On exit the RNG state in \code{.Random.seed} is reset to the value that applied when the function was called.
}

For NULL seed input, the saved RNGstate may be used to reset .Random.seed (see Examples).

}

\section{Use of random numbers in \pkg{secr}}{

Many functions in \pkg{secr} call on random numbers, sometimes in unexpected places. For example \code{\link{autoini}} selects a random sample to thin points and speed computation. In most functions there is no provision for direct control of the random number state: users won't usually care, and if they do then \code{set.seed} may be called for the particular \R session. 

However, control of the RNG seed is required for reproducible data generation in simulation functions. These functions typically have a `seed' argument that is used internally in a call to \code{set.seed}. Handling of seeds in the simulation functions of \pkg{secr} largely follows \code{stats::simulate} as described in the preceding section.

The relevant functions are --
\tabular{llll}{
Function \tab Default \tab Saved attribute \tab Note \cr
\code{\link{randomHabitat}} \tab NULL \tab seed or RNGstate\cr
\code{\link{secr.test}} \tab NULL \tab seed or RNGstate \tab calls and retains seed from \code{simulate.secr} \cr
\code{\link{sim.capthist}} \tab NULL \tab seed or RNGstate \cr
\code{\link{sim.resight}} \tab NULL \tab seed or RNGstate \tab Seed may be passed in \dots argument \cr
\code{\link{sim.popn}} \tab NULL \tab seed or RNGstate \cr
\code{\link{sim.secr}} \tab NULL \tab seed or RNGstate \cr
\code{\link{simulate.secr}} \tab NULL \tab seed or RNGstate \tab S3 method called by \code{sim.secr}\cr
}

Setting \code{seed = NULL} in any of these functions has the effect of continuing the existing random number stream; it is not the same as calling \code{set.seed(NULL)}. \code{\link{par.secr.fit}} does not follow this model.

}

\section{Parallel processing}{

Two models are used for parallel processing in \pkg{secr}, corresponding to multi-threading with package \pkg{RcppParallel} (e.g. \code{\link{secr.fit}}) and parallel cores in package \pkg{parallel} (e.g. \code{\link{chat.nk}}).

In the \pkg{parallel} model the L'Ecuyer pseudorandom generator is used to provide a separate random number stream for each core (see \code{\link{clusterSetRNGStream}}).

When using \pkg{Rcpp} the state of the random number generator is set in C++ with the call 
\code{RNGScope scope;}that automatically resets the state of the generator on exit (Eddelbuettel 2013 p. 115).

Random number streams in separate \pkg{RcppParallel} threads are (probably) not independent. Thus there are potential issues with RNG calls in multi-threaded code. However, in \pkg{secr} 4.0 all RNG calls in C++ code are outside multi-threaded contexts, with the exception of simulations allowing for overdispersion in mark--resight estimates (Rcpp exported function sightingchatcpp). The implications for mark-resight estimates have not been explored, and it is unclear whether more elaborate solutions are needed.

}


\seealso{

  \code{\link{Parallel}},
  \code{\link{set.seed}}, 
  \code{\link{simulate}}, 
  \code{\link{sim.capthist}}, 
  \code{\link{sim.popn}}, 
  \code{\link{sim.resight}},
  \code{\link{secr.test}},
  \code{\link{simulate.secr}}

}

\references{

Eddelbuettel, D. 2013. Seamless R and C++ integration with Rcpp. Springer.

}

\examples{

\dontrun{

lmfit <- lm(speed ~ dist, data = cars)

## 1. NULL seed
r1 <- simulate(lmfit, seed = NULL)
r2 <- simulate(lmfit, seed = NULL)
## restore RNGstate, assuming RNGkind unchanged
.Random.seed <- attr(r1, "seed")
r3 <- simulate(lmfit, seed = NULL)
r1[1:6,1]
r2[1:6,1]
r3[1:6,1]

## 2. explicit seed
r4 <- simulate(lmfit, seed = 123)
r5 <- simulate(lmfit, seed = attr(r4, "seed"))
r4[1:6,1]
r5[1:6,1]

}

}

\keyword{ datagen }
