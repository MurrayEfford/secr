\name{occasionKey}
\alias{occasionKey}

\title{Key to Petal Plot}

\description{Displays a graphic key to the occasions corresponding to
  petals in a petal plot.}

\usage{
occasionKey(capthist, noccasions, rad = 3, x, y, px = 0.9, py = 0.9, 
                         title = 'Occasion', ...)
}

\arguments{
  \item{capthist}{single-session capthist object}
  \item{noccasions}{number of petals (if \code{capthist} not provided)}
  \item{rad}{distance of petal centre from key centre}
  \item{x}{numeric x coordinate for centre of key}
  \item{y}{numeric y coordinate for centre of key}
  \item{px}{x position as fraction of user coordinates}
  \item{py}{y position as fraction of user coordinates}
  \item{title}{character}
  \item{\dots}{other arguments passed to \code{\link{text}}}
}

\details{
  
Either \code{capthist} or \code{noccasions} is required. It is assumed that a plot exists.

Graphic arguments in \dots are applied to both the title and the occasion numbers.

}

\value{

The key will be added to an existing plot. No value is returned.

}

\seealso{
  
  \code{\link{plot.capthist}}
  
}

\examples{

plot(captdata, border = 50)
occasionKey(captdata, rad = 8, cex = 0.8)

}

\keyword{ hplot }
