\name{utility}
\alias{getMeanSD}
\alias{maskarea}
\alias{masklength}
\alias{edist}
\alias{nedist}

\title{Utility Functions}

\description{
  
Minor functions.

}

\usage{
getMeanSD(xy)
maskarea(mask, sessnum = 1)
masklength(mask, sessnum = 1)
edist(xy1, xy2)
nedist(xy1, xy2, mask, inf = Inf, ...)
rlnormCV(n, mean, CV)
}

\arguments{ 
  \item{xy}{2-column matrix or dataframe}
  \item{xy1}{2-column matrix or dataframe}
  \item{xy2}{2-column matrix or dataframe}
  \item{mask}{mask or linearmask object}
  \item{sessnum}{integer; for multi-session masks, the number of the session}
  \item{inf}{numeric value to use for +infinity}
  \item{\dots}{other arguments for \code{\link[gdistance]{transition}}}
  \item{n}{number of observations}
  \item{mean}{mean on natural scale}
  \item{CV}{coefficient of variation on natural scale}
}

\details{

  \code{getmeanSD} is used by \code{\link{make.mask}} to standardize
  mask coordinates.
  
  For \code{masklength} the input should be a linear mask from \pkg{secrlinear}.

  \code{edist} computes the Euclidean distance between each point in xy1
  and each point in xy2. (This duplicates the functionality of `rdist'
  in package \pkg{fields}).

  \code{nedist} computes the non-Euclidean distance between each point
  in xy1 and each point in xy2, in two dimensions. The calculation uses
  \pkg{gdistance} (van Etten 2017; see also Csardi & Nepusz 2006): a
  transition layer is formed representing the connections between
  adjacent points in \code{mask}. By default, points within a 16-point
  neighbourhood are considered `adjacent'. Distances are obtained by
  Dijkstra's (1959) algorithm as least cost paths through the graph of
  all points in the mask.

  \code{nedist} has some subtle options. If `mask' is missing then the
  transition layer will be formed from `xy2'. If `mask' has a covariate
  named `noneuc' then this will be used to weight distances. The \dots
  argument of \code{nedist} allows the user to vary arguments of
  \code{\link[gdistance]{transition}} (defaults transitionFunction =
  mean and directions = 16). Be warned this can lead to unexpected
  results! Point pairs that are completely separated receive the
  distance +Inf unless a finite value is provided for the argument
  `inf'. See
  \href{https://www.otago.ac.nz/density/pdfs/secr-noneuclidean.pdf}{secr-noneuclidean.pdf}
  for uses of \code{nedist}.
  
  \code{rlnormCV} is a wrapper for \code{\link{rlnorm}} that computes its 
  meanlog and sdlog arguments from the mean and CV on the natural scale:
  
  \eqn{\mbox{sdlog} = \sqrt{\log(\mbox{CV}^2 + 1)}}, and  
  
  \eqn{\mbox{meanlog} = \log(\mbox{mean}) - \mbox{sdlog}^2/2}.
  
}

\value{

  For \code{getMeanSD}, a dataframe with columns `x' and `y' and two
  rows, mean and SD.

  For \code{maskarea}, the summed area of mask cells in hectares (ha).

  For \code{masklength}, the summed length of mask cells in kilometers (km).
  
  For \code{edist} and \code{nedist}, a matrix with dim = c(nrow(xy1), nrow(xy2)).

  For \code{rlnormCV} a vector of random deviates.
}

\references{

Dijkstra, E. W. (1959) A note on two problems in connexion with
graphs. \emph{Numerische Mathematik}, \bold{1}, 269--271.

Csardi, G. and Nepusz, T. (2006) The igraph software package for complex
network research. \emph{InterJournal}, \bold{1695}. \url{https://igraph.org}

van Etten, J. (2017) R package gdistance: Distances and routes on
  geographical grids. \emph{Journal of Statistical Software}, \bold{76(1)}, 1--21.
  \doi{10.18637/jss.v076.i13}

}

\examples{

getMeanSD(possummask)

}

