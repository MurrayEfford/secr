\name{discretize}
\alias{discretize}

\title{Rasterize Area Search or Transect Data}

\description{

It is sometimes useful to re-cast area-search (polygon or polygonX) data as if it was from a
set of closely spaced point detectors, i.e. to rasterize the detection locations. This function makes that conversion. Each polygon detector in the input is replaced by a number of point detectors, each representing a square pixel. Detections are mapped to the new detectors on the basis of their x-y coordinates.
  
If \code{object} contains transect data the problem is passed to \code{\link{snip}} and \code{reduce.capthist}.
  
}

\usage{

discretize(object, spacing = 5, outputdetector = c("proximity", "count", "multi"), 
    tol = 0.001, cell.overlap = FALSE, type = c("centre","any", "all"), ...)

}

\arguments{
  \item{object}{secr capthist or traps object}
  \item{spacing}{numeric spacing between point detectors in metres }
  \item{outputdetector}{character output detector type}
  \item{tol}{numeric fractional inflation of perimeter (see Details)}
  \item{cell.overlap}{logical; if TRUE the area of overlap is stored in usage attribute }
  \item{type}{character; see Details}
  \item{\dots}{other arguments passed to \code{\link{snip}} if object is transect}
}

\details{
  
  The input should have detector type `polygon' or `polygonX'.

  A new array of equally spaced detectors is generated within
  each polygon of the input, inflated radially by 1 + tol to avoid some
  inclusion problems. The origin of the superimposed grid is fixed automatically. If \code{type = "centre"} detectors are included if they lie within the (inflated) polygon. Otherwise, the decision on whether to include a candidate new detector is based on the corner vertices of the cell around the detector (side = \code{spacing}); \code{type = "any"} and  \code{type = "all"} have the obvious meanings. 

  \code{tol} may be negative, in which case the array(s) will be shrunk relative
  to the polygon(s). 
  
  For irregular polygons the edge cells in the output may be only partially contained within the polygon they represent. Set \code{cell.overlap = TRUE} to retain the proportion of overlap as the `usage' of the new traps object. This can take a few seconds to compute. If `usage' is already defined then the new `usage' is the old multiplied by the proportion of overlap.
  
  Combining \code{type = "any"} and \code{cell.overlap = TRUE} with \code{tol} > 0 can have the odd effect of including some marginal detectors that are assigned zero usage.
  
  With \code{type = "any"}, the sum of the overlap proportions times cell area is equal to the area of the polygons.

}

\value{
  
  A capthist or traps object of the requested detector type, but otherwise
  carrying forward all attributes of the input. The embedded traps object has a 
  factor covariate `polyID' recording the polygon to which each point detector relates.
  
}

\note{
  
  Consider the likely number of detectors in the output before you start. 
  
}

\seealso{
  
  \code{\link{reduce.capthist}}, \code{\link{snip}}
  
}

\examples{

\dontrun{

## generate some polygon data
pol <- make.poly()
CH <- sim.capthist(pol, popn = list(D = 30), detectfn = 'HHN', 
   detectpar = list(lambda0 = 0.3))
plot(CH, border = 10, gridl = FALSE, varycol = FALSE)

## discretize and plot
CH1 <- discretize(CH, spacing = 10, output = 'count')
plot(CH1, add = TRUE, cappar = list(col = 'orange'), varycol =
    FALSE, rad = 0)
plot(traps(CH1), add = TRUE)
# overlay cell boundaries
plot(as.mask(traps(CH1)), dots = FALSE, col = NA, meshcol = 'green', 
    add = TRUE)

## show how detections are snapped to new detectors
newxy <- traps(CH1)[nearesttrap(xy(CH),traps(CH1)),]
segments(xy(CH)[,1], xy(CH)[,2], newxy[,1], newxy[,2])

plot(traps(CH), add = TRUE)  # original polygon

## Incomplete overlap

pol <- rotate(make.poly(), 45)
CH2 <- sim.capthist(pol, popn = list(D = 30), detectfn = 'HHN', 
   detectpar = list(lambda0 = 0.3))
plot(CH2, border = 10, gridl = FALSE, varycol = FALSE)
CH3 <- discretize(CH2, spacing = 10, output = 'count', type = 'any', 
    cell.overlap = TRUE, tol=0.05)
    
plot(CH3, add = TRUE, cappar = list(col = 'orange'), varycol =
    FALSE, rad = 0)
plot(traps(CH3), add = TRUE)

# overlay cell boundaries and usage
msk <- as.mask(traps(CH3))
covariates(msk) <- data.frame(usage = usage(traps(CH3))[,1])
plot(msk, dots = FALSE, cov='usage', meshcol = 'green', 
    add = TRUE)
    
}

}

\keyword{ manip }

