\name{plot.popn}
\alias{plot.popn}

\title{ Plot Population Object }

\description{

Display animal locations from a \code{popn} object.

}

\usage{
\method{plot}{popn}(x, add = FALSE, frame = TRUE,  
    circles = NULL, collapse = FALSE, seqcol = NULL, ...)
}

\arguments{
  \item{x}{ object of class \code{popn} }
  \item{add}{ logical to add points to an existing plot }
  \item{frame}{ logical to add frame or polygon within which points were
      simulated}
  \item{circles}{ vector giving the radii if circles are to be plotted }
  \item{collapse}{ logical; if TRUE then multiple sessions are overlaid }
  \item{seqcol}{ color used for first detection when collapse = TRUE (optional) }
  \item{\dots}{ arguments passed to \code{eqscplot} and \code{points} or
  \code{symbols}}
  
}

\details{
  If \code{circles} is provided then a circle of the given radius is
  plotted for each animal using the \code{symbols} function. The
  arguments \code{fg} and \code{bg} may be used to control the colour of
  the perimeter and the fill of each circle (see Examples).
  
  For a multi-session popn with \code{\link{turnover}}, \code{collapse = TRUE} 
  allows successive locations to be joined with (\code{type = 'o'} or \code{type = 'l'}).
  
  \code{seqcol} may be a single color, a vector of colours (one per session), 
  or a vector of two colours, one for the first and one for all later sessions 
  in which each animal was detected. 
  
  If \dots includes 'col' then 'collapse' must be specified to avoid confusion, 
  even for single-session data (see Examples).
  
}

\seealso{ \code{\link{popn}}, \code{\link{sim.popn}} }

\examples{

temppopn <- sim.popn(D = 5, expand.grid(
    x = c(0,100), y = c(0,100)))
# specify collapse to avoid partial match of col    
plot(temppopn, pch = 16, collapse = FALSE, col = "blue")

plot(temppopn, circles = 20, bg = "tan", fg = "white")
plot(temppopn, pch = 16, cex = 0.5, add = TRUE)

}

