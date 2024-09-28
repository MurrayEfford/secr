\name{popn}
\alias{popn}

\title{ Population Object }

\description{
Encapsulate the locations of a set of individual animals.
}

\details{
  
An object of class \code{popn} records the locations of a set of
individuals, together with ancillary data such as their sex.  Often used
for a realisation of a spatial point process (e.g. homogeneous Poisson)
with known density (intensity).  Locations are stored in a data frame
with columns `x' and `y'.

A \code{popn} object has attributes
\tabular{ll}{
covariates \tab data frame with numeric, factor or character variables to be used as individual covariates \cr
model2D    \tab 2-D distribution ("poisson", "cluster", "IHP", "linear" etc.) \cr
Ndist      \tab distribution of number of individuals ("poisson", "fixed") \cr
boundingbox \tab data frame of 4 rows, the vertices of the rectangular area \cr
}

The number of rows in \code{covariates} must match the length of
\code{x} and \code{y}. See \code{\link{sim.popn}} for more information
on \code{Ndist} and \code{model2D}.

}

\note{ The \code{popn} class is used only occasionally: it is not central to spatially explicit capture recapture. }

\seealso{ \code{\link{sim.popn}}, \code{\link{plot.popn}}, \code{\link{transformations}} }


