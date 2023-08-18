\name{as.popn}
\alias{as.popn}

\title{Coerce ppp object to popn}

\description{

This function converts a \pkg{spatstat} "ppp" object (Baddeley et al. 2015), making it easier to use the simulation capability of \pkg{spatstat} in \pkg{secr}.

}

\usage{

as.popn(x)

}

\arguments{

  \item{x}{an object of class 'ppp'}

}

\details{

  Not all attributes are carried over.
  
}

\value{

An object of class c("popn", "data.frame") with attribute "boundingbox". The attribute "Lambda"  (\pkg{spatstat} class "im") is also carried over if present (used for the intensity surface of LGCP simulations).

}

\seealso{

  \code{\link{sim.popn}}, \code{\link{popn}}

}

\references{

Baddeley, A., Rubak, E., and Turner, R. 2015. Spatial Point Patterns: Methodology and Applications with R. Chapman and Hall/CRC Press, London. ISBN 9781482210200, https://www.routledge.com/Spatial-Point-Patterns-Methodology-and-Applications-with-R/Baddeley-Rubak-Turner/p/book/9781482210200/. 

}

\keyword{manip}
