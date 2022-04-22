\name{plotMaskEdge}
\alias{plotMaskEdge}

\title{Outline Around Mask Cells}

\description{
  
Plots the outer edge of a  mask.

}
\usage{

plotMaskEdge(mask, plt = TRUE, add = FALSE, ...)

}

\arguments{
  \item{mask}{secr habitat mask object}
  \item{plt}{logical; if TRUE the edge is plotted}
  \item{add}{logical; if TRUE the line is added to an existing plot}
  \item{\dots}{other line plotting arguments passed to \code{\link{segments}}}
}

\details{
  
  May be slow. 
  
}

\value{
  A numeric matrix of 4 columns is returned invisibly. The columns may be used as arguments x0, y0, x1, y1 in a call to segments().
}

\note{

A bug in \pkg{secr} <3.2.2 caused some internal lines to appear when the mask spacing was not an integer.

}

\examples{

\dontrun{
plot(possummask)
plotMaskEdge (possummask, add = TRUE)
}

}

\keyword{ hplot }
