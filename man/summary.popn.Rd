\name{summary.popn}
\alias{summary.popn}
\alias{print.summary.popn}

\title{ Summarise Simulated Population }

\description{

 Concise summary of a \code{popn} object.

}

\usage{
\method{summary}{popn}(object, collapse = FALSE, ...)
\method{print}{summary.popn}(x, ...)

}

\arguments{
    \item{object}{ \code{popn} object }
    \item{collapse}{logical; if TRUE multi-session popn objects are treated as a single open population}
  \item{x}{ \code{summary.popn} object }
  \item{\dots}{ other arguments (not used) }
}

\details{

By default each component of a multisession object is summarised separately. If \code{collapse = TRUE} then turnover and movements are collated across sessions, matching individuals by rownames.

}

\value{

For \code{summary.popn}, an object of class `summary.popn' with various components. For a multisession object and collapse = TRUE the descriptors include the numbers of new individuals (recruits) and lost individuals (deaths), and matrices showing the status of each animal in each session (`status' codes 0 not recruited yet; 1 alive; -1 dead) and movement from previous session if alive then (`movements').

}

\seealso{

\code{\link{sim.popn}}

}

\examples{


grid <- make.grid(8,8)
turnover <- list(phi = 0.8, lambda = 1)
pop <- sim.popn(Nbuffer = 200, core = grid, buffer = 200, Ndist = 'fixed', 
    nsessions = 5, details = turnover)
summary(pop, collapse = TRUE)

}

\keyword{ datagen }