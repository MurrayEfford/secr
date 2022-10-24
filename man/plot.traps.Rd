\name{plot.traps}
\alias{plot.traps}

\title{ Plot traps Object }

\description{
  Map the locations of detectors (traps). 
}

\usage{
\method{plot}{traps}(x, border = 100, label = FALSE, offset = c(6,6), add = FALSE,
  hidetr = FALSE, detpar = list(), txtpar = list(), bg = "white",
  gridlines = !add, gridspace = 100, gridcol = "grey", 
  markused = FALSE, markvarying = FALSE, markvertices = FALSE,
  labelclusters = FALSE, frame = NULL, ...)
}

\arguments{
  \item{x}{ a \code{traps} object }
  \item{border}{ width of blank margin around the outermost detectors }
  \item{label}{ logical indicating whether a text label should appear by each detector }
  \item{offset}{ vector displacement of label from point on x and y axes }
  \item{add}{ logical to add detectors to an existing plot } 
  \item{hidetr}{ logical to suppress plotting of detectors }
  \item{detpar}{ list of named graphical parameters for detectors (passed to \code{par}) }
  \item{txtpar}{ list of named graphical parameters for labels (passed to \code{par}) }
  \item{bg}{ background colour }
  \item{gridlines}{ logical for plotting grid lines }
  \item{gridspace}{ spacing of gridlines }
  \item{gridcol}{ colour of gridlines }
  \item{markused}{ logical to distinguish detectors used on at least one
  occasion}
  \item{markvarying}{ logical to distinguish detectors whose usage varies among occasions }
  \item{markvertices}{ logical or 0,1,2 for plotting transect or polygon points }
  \item{labelclusters}{ logical to label clusters}
  \item{frame}{data defining a boundary polygon (see \code{\link{boundarytoSF}}  }
  \item{\dots}{ arguments to pass to \code{eqscplot} }
}
\details{
  
\code{offset} may also be a scalar value for equal displacement on the x
and y axes. The \code{hidetr} option is most likely to be used when
plot.traps is called by plot.capthist. See \code{\link{par}} and
\code{\link{colours}} for more information on setting graphical
parameters. The initial values of graphical parameters are restored on
exit.

Axes are not labeled. Use \code{\link{axis}} and \code{\link{mtext}} if
necessary.

\code{markvertices} determines whether the vertices of each transect or
polygon will be emphasised by overplotting a point symbol
(detpar$pch). Value may be logical (TRUE, FALSE) or integer (0 = no
points, 1 = terminal vertices only, 2 = all vertices).

From 4.4.0, polygon detectors are shaded with detpar$col and outlined (border) with detpar$fg. Use detpar$col = NA for no shading (transparent polygons).

\code{labelclusters} requires \code{x} to have attributes `clusterID' and
`clustertrap'.

A boundary polygon is plotted in black if \code{frame} is not NULL.

}

\value{
None
}

\seealso{ \code{\link{plot}}, \code{\link{traps}},
  \code{\link{clusterID}}  }

\examples{
 temptrap <- make.grid()
 plot (temptrap, detpar = list(pch = 16, col = "blue"), 
    label = TRUE, offset = 7)
}

\keyword{ hplot }

