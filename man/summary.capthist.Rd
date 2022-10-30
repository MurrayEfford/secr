\name{summary.capthist}
\alias{summary.capthist}
\alias{print.summary.capthist}
\alias{counts}

\title{ Summarise Detections }

\description{
  Concise description of \code{capthist} object.
}
\usage{
\method{summary}{capthist}(object, terse = FALSE, moves = FALSE, tpa = FALSE, ...)

\method{print}{summary.capthist}(x, ...)

counts(CHlist, counts = "M(t+1)")

}

\arguments{
  \item{object}{ \code{\link{capthist}} object }
  \item{terse}{ logical; if TRUE return only summary counts }
  \item{moves}{ logical; if TRUE then summary includes detected movements }
  \item{tpa}{ logical; if TRUE then summary includes number of detectors per animal}
  \item{x}{ \code{summary.capthist} object }
  \item{\dots}{ arguments passed to other functions }

  \item{CHlist}{ capthist object, especially a multi-session object }
  \item{counts}{ character vector of count names }

}

\value{

From \code{summary.capthist}, an object of class
\code{summary.capthist}, a list with at least these components

  \item{detector }{\code{\link{detector}} type ("single", "multi", "proximity" etc.)}
  \item{ndetector }{number of detectors }
  \item{xrange }{range of x coordinates of detectors }
  \item{yrange }{range of y coordinates of detectors }
  \item{spacing }{mean distance from each trap to nearest other trap }
  \item{counts }{matrix of summary counts (rows) by occasion (columns).  See Details.}
  \item{dbar }{mean recapture distance }
  \item{RPSV }{root pooled spatial variance }

or, when \code{terse = TRUE}, a vector (single session) or dataframe (multiple sessions) of counts (Occasions, Detections, Animals, Detectors, and optionally Moves and Animals2).

A summary of individual covariates is provided if these are present (from \pkg{secr} 4.0.1).

A summary of interference/non-target captures is provided if there is a \link{nontarget} attribute (from \pkg{secr} 4.5.5).
  
}

\details{
These counts are reported by \code{summary.capthist}  
\tabular{ll}{
  n \tab number of individuals detected on each occasion \cr
  u \tab number of individuals detected for the first time on each occasion \cr
  f \tab number of individuals detected exactly f times \cr
  M(t+1) \tab cumulative number of individuals detected \cr
  losses \tab number of individuals reported as not released on each occasion \cr
  detections \tab number of detections, including within-occasion `recaptures' \cr
  traps visited \tab number of detectors at which at least one detection was recorded \cr
  traps set \tab number of detectors, excluding any `not set' in \code{usage} attribute of \code{traps} attribute \cr
}
The last two rows are dropped if the data are nonspatial (object has no traps attribute).

Movements are as reported by \code{\link{moves}}. When \code{terse = TRUE} the number of non-zero moves is reported. The temporal sequence of detections at `proximity' and `count' detectors is not recorded in the capthist object, so the movement statistics are not to be taken too seriously. The problem is minimised when detections are sparse (seldom more than one per animal per occasion), and does not occur with `single` or `multi` detectors.

The `tpa' option provides the frequency distribution of detectors per animal. When \code{terse = TRUE} the number of animals at >= 2 detectors is reported ('Animal2').

\code{counts} may be used to return the specified counts in a compact 
session x occasion table. If more than one count is named then a list is returned
with one component for each type of count.

}

\seealso{ \code{\link{dbar}}, \code{\link{RPSV}}, \code{\link{capthist}} }

\examples{
temptrap <- make.grid(nx = 5, ny = 3)
summary(sim.capthist(temptrap))
summary(sim.capthist(temptrap))$counts["n",]
summary(captdata, moves = TRUE)
}

\keyword{ models }

