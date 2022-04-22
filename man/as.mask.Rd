\name{as.mask}
\alias{as.mask}

\title{Coerce traps object to mask}

\description{This function is used primarily for plotting covariates, for which the plot.mask function has greater functionality than \code{plot.traps}.  It also generates pretty maps of grid cells. }

\usage{

as.mask(x)

}

\arguments{
  \item{x}{an object of class 'traps'}
}

\details{

  A mask derived by coercion with \code{as.mask} may behave
  unpredictably e.g., in \code{secr.fit}.
  
}

\value{

  If \code{x} is a single-session traps object --
  
  an object of class c("mask", "data.frame")

  If \code{x} is a multi-session traps object --
  
  an object of class c("mask", "list"), for which each component is a
  single-session mask.
  
}

\seealso{

  \code{\link{make.mask}},  \code{\link{plot.mask}},  
  \link{mask}, \link{traps}

}

\examples{

plot(as.mask(traps(captdata)), dots = FALSE, meshcol = "black")
plot(traps(captdata), add = TRUE)

}

\keyword{manip}
