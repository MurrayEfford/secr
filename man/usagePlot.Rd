\name{usagePlot}
\alias{usagePlot}
\alias{sightingPlot}

\title{
Plot usage, detections or sightings.
}

\description{

\code{usagePlot} displays variation in effort (usage) over detectors as a
bubble plot (circles with radius scaled so that area is proportional to
effort).

\code{sightingPlot} displays spatial variation in the number of sightings at each detector as a
bubble plot (circles with radius scaled so that area is proportional to either the average number per occasion or the total over occasions.

}

\usage{

usagePlot(object, add = FALSE, occasions = NULL, col = "black", fill =
FALSE, scale = 2, metres = TRUE, rad = 5, ...)

sightingPlot(object, type = c("Detections", "Tu", "Tm", "Tn"), add = FALSE, 
occasions = "ALL", mean = TRUE, col = "black", fill = FALSE, scale = 2, 
metres = TRUE, dropunused = TRUE, title = type, legend = c(1, 2, 4, 8), 
px = 0.95, py = 0.95, ...)

}

\arguments{
  \item{object}{traps object with usage attribute}
  \item{add}{logical; if FALSE plot.traps is called to create a base plot}
  \item{occasions}{integer number(s) of the occasion(s) for which effort is
  plotted, "ALL", or NULL}
  \item{col}{character or integer colour value}
  \item{fill}{logical; if TRUE the circle is filled with the line colour}
  \item{scale}{numeric value used to scale radius}
  \item{metres}{logical; if TRUE scale is a value in metres (see Details)}
  \item{rad}{numeric; radial displacement of symbol centre for each occasion from
  true detector location (metres) }
  \item{\dots}{other arguments passed to plot.traps}
  
  \item{type}{character to choose among sighting types and detections of marked animals}
  \item{mean}{logical; if TRUE then the plotted value is the average over occasions, otherwise the sum}
  \item{dropunused}{logical; if TRUE then detectors are omitted when they were unused on \code{occasions}}
  \item{title}{character}
  \item{legend}{numeric values for which legend circles will be drawn}
  \item{px}{legend x position as fraction of user coordinates}
  \item{py}{legend y position as fraction of user coordinates}
}

\details{

The behaviour of \code{usagePlot} is described first. By default (\code{occasion = NULL}) circles representing usage on each
  occasion are plotted around the detector location at distance
  \code{rad}, as in the petal plot of
  \code{\link{plot.capthist}}. Otherwise, the usage on a single
  specified occasion, or summed over occasions
  (\code{length(occasion)>1}, or \code{occasion = "ALL"}), is plotted as
  a circle centred at the detector location.

Package \pkg{sp} provides an alternative to \code{usagePlot} (see Examples).

\code{sightingPlot} may be used to display either detections of marked animals (whether or not \code{occasions} refers to sighting occasions) or any of the sighting attributes (unmarked sightings `Tu', marked, unidentified sightings `Tm', or other uncertain sightings `Tn').

If \code{py} is of length 2 then the values determine the vertical spread of symbols in the legend.

For both functions --

The \code{metres} argument switches between two methods. If \code{metres
= TRUE}, the \code{symbols} function is used with \code{inches = FALSE} to plot
circles with radius scaled in the units of \code{object} (i.e. metres;
\code{scale} is then the radius in metres of the symbol for a detector
with usage = 1.0).  Otherwise, plotting uses \code{points}; this has the
advantage of producing better filled circles, but a suitable value of
scale must be found by trial and error.
  
}

\value{

  No value is returned by \code{usagePlot}.
  
  \code{sightingPlot} invisibly returns a `traps' object with a covariate `f' holding the plotted values.
  
}

\seealso{
  
  \code{\link{usage}}, \code{\link{symbols}}, \code{\link[sp]{bubble}},
  \link{sightings}
  
}

\examples{

simgrid <- make.grid(nx = 10, ny = 10, detector = "proximity")
usage(simgrid) <- matrix(rep(1:10, 50), nrow = 100, ncol = 5)
usagePlot(simgrid, border = 20, scale = 1.5, fill = FALSE,
  metres = FALSE)

# It is hard to get the legend just right
# here is one attempt
legend (x = -50, y = 185, legend = c(1,2,5,10), pch = 1, pt.cex =
  c(1,2,5,10)^0.5 * 1.5, x.intersp = 3, y.intersp = 1.8, adj = 1,
  bty = "n", title = "Usage")

usagePlot(simgrid, occasion = NULL, border = 20, scale = 1.5, fill = FALSE,
  metres = FALSE)

\dontrun{
# bubble plot in package 'sp'
library(sp)
simgrid$usage <- usage(simgrid)[,1]  ## occasion 1
class(simgrid) <- "data.frame"
coordinates(simgrid) <- c("x","y")
bubble(simgrid)   
}

}


\keyword{ hplot }
