\name{deleteMaskPoints}
\alias{deleteMaskPoints}

\title{Edit Mask Points}

\description{Mask points may be removed by one of three methods:
  clicking on points, clicking on vertices to define a polygon from
  which points will be removed, or specifying a polygon to which the
  mask will be clipped.
  }

\usage{

deleteMaskPoints(mask, onebyone = TRUE, add = FALSE, poly = NULL,
 poly.habitat = FALSE, ...)

}

\arguments{
  \item{mask}{secr mask object}
  \item{onebyone}{logical; see Details}
  \item{add}{logical; if true then the initial mask plot will be added
  to an existing plot}
  \item{poly}{polygon defining habitat or non-habitat as described in
    \code{\link{make.mask}}}
  \item{poly.habitat}{logical; if TRUE polygon represents habitat}
  \item{\dots}{other arguments to plot.mask}
}

\details{

  The default method (onebyone = TRUE, poly = NULL) is to click on each
  point to be removed. The nearest mask point will be selected.

  Setting onebyone = FALSE allows the user to click on the vertices of a
  polygon within which all points are to be removed (the default) or
  retained (\code{poly.habitat = TRUE}). Vertices need not
  coincide with mask points.

  Defining \code{poly} here is equivalent to calling \code{make.mask}
  with \code{poly} defined. \code{poly} one of the several types described 
  in \code{\link{boundarytoSF}}. Whether \code{poly} represents habitat or
  non-habitat is toggled with \code{poly.habitat} -- the default here
  differs from \code{\link{make.mask}}.

}

\value{
  
  A mask object, usually with fewer points than the input mask.
  
}

\seealso{
  
  \code{\link{make.mask}}, \code{\link{subset.mask}}
  
}

\examples{

if (interactive()) {
    mask0 <- make.mask (traps(captdata))
    ## Method 1 - click on each point to remove
    mask1 <- deleteMaskPoints (mask0)
    ## Method 2 - click on vertices of removal polygon
    mask2 <- deleteMaskPoints (mask0, onebyone = FALSE)
    ## Method 3 - predefined removal polygon
    plot(captdata)
    poly1 <- locator(5)
    mask3 <- deleteMaskPoints (mask0, poly = poly1)
}

}
\keyword{ manip }
