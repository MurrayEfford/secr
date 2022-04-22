\name{plot.capthist}
\alias{plot.capthist}
\alias{plotMCP}

\title{ Plot Detection Histories }

\description{
  Display a plot of detection (capture) histories or telemetry data over a map of the detectors.
}

\usage{
\method{plot}{capthist}(x, rad = 5, hidetraps = FALSE, tracks = FALSE,
    title = TRUE, subtitle = TRUE, add = FALSE, varycol = TRUE, 
    icolours = NULL, randcol = FALSE, lab1cap = FALSE, laboffset = 4, ncap = FALSE, 
    splitocc = NULL, col2 = "green", type = c("petal", "n.per.detector", "n.per.cluster",
        "sightings", "centres", "telemetry"),
    cappar = list(cex = 1.3, pch = 16, col = "blue"), 
    trkpar = list(col = "blue", lwd = 1), 
    labpar = list(cex = 0.7, col = "black"), ...)

plotMCP(x, add = FALSE, col = "black", fill = NA, lab1cap = FALSE, laboffset = 4,
    ncap = FALSE, ...)

}

\arguments{
  \item{x}{ an object of class \code{capthist} }
  \item{rad}{ radial displacement of dot indicating each capture event from the detector location (used to separate overlapping points) }
  \item{hidetraps}{ logical indicating whether trap locations should be displayed }
  \item{tracks}{ logical indicating whether consecutive locations of individual animals should be joined by a line }
  \item{title}{ logical or character string for title }
  \item{subtitle}{ logical or character string for subtitle }
  \item{add}{ logical for whether to add to existing plot }
  \item{varycol}{ logical for whether to distinguish individuals by colour }
  \item{icolours}{ vector of individual colours (when \code{varycol} =
  TRUE), or colour scale (non-petal plots) }
  \item{randcol}{ logical to use random colours (\code{varycol} = TRUE) }
  \item{lab1cap}{ logical for whether to label the first capture of each animal }
  \item{laboffset}{ distance by which to offset labels from points }
  \item{ncap}{ logical to display the number of detections per trap per occasion }
  \item{splitocc}{ optional occasion from which second colour is to be used  }
  \item{col2}{ second colour (used with \code{splitocc}) }
  \item{type}{ character string ("petal", "n.per.detector" or "n.per.cluster") }
  \item{cappar}{ list of named graphical parameters for detections (passed to \code{par}) }
  \item{trkpar}{ list of named graphical parameters for tracks (passed to \code{par}) }
  \item{labpar}{ list of named graphical parameters for labels (passed to \code{par}) }
  \item{\dots}{ arguments to be passed to \code{plot.traps} }
  \item{col}{vector of line colour numbers or names (plotMCP only)}
  \item{fill}{vector of fill colour numbers or names (plotMCP only)}

}

\details{ 

By default, a `petal' plot is generated in the style of Density (Efford 2012)
using \code{eqscplot} from the MASS library. 

If \code{type =
"n.per.detector"} or \code{type = "n.per.cluster"} the result is a
colour-coded plot of the number of individuals at each unit, pooled over
occasions. 

If \code{type = "sightings"} the sightings of unmarked animals are 
displayed on a petal-like plot (requires mark-resight data) (see also \code{\link{sightingPlot}}).

If \code{type = "centres"} then a single point is plotted for each animal, jittered on each axis by a random amount (limits +/- \code{rad}/2).

If \code{type = "telemetry"} and the `telemetryxy' attribute is not NULL then the telemetry locations are plotted.

If \code{title} = FALSE no title is displayed; if \code{title} = TRUE,
the session identifer is used for the title.

If \code{subtitle} = FALSE no subtitle is displayed; if \code{subtitle}
= TRUE, the subtitle gives the numbers of occasions, detections and
individuals.

If \code{x} is a multi-session capthist object then a separate plot is
produced for each session. Use \code{par(mfrow = c(nr, nc))} to allow a
grid of plots to be displayed simultaneously (nr rows x nc columns).

These arguments are used only for petal plots: \code{rad},
\code{tracks}, \code{varycol}, \code{randcol}, \code{lab1cap},
\code{laboffset}, \code{ncap}, \code{splitocc}, \code{col2},
\code{trkpar}, and \code{labpar}. Call \code{\link{occasionKey}} 
to add a key to the petals.

If \code{icolours = NULL} and \code{varycol = TRUE} then a vector of
colours is generated automatically as topo.colors((nrow(x)+1) * 1.5).
If there are too few values in \code{icolours} for the number of
individuals then colours will be re-used. 

\code{plotMCP} plots minimum convex polygons of individual location
data over a base plot of detector locations. Usually the data are
telemetry locations in the xylist attribute of the capthist
object; if this is not present and \code{x} is a polygon search
capthist then the individual xy data are plotted.

}

\value{
  
  For \code{type = "petal"}, the number of detections in \code{x}.
  For \code{type = "sightings"}, the number of sightings of unmarked animals in \code{x}.
  For \code{type = "n.per.detector"} or \code{type = "n.per.cluster"}, a
  dataframe with data for a legend (see Examples).

  \code{plotMCP} invisibly returns a list in which each component is a
  2-column (x,y) dataframe of boundary coordinates for one individual.

}

\references{ 

Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially explicit
capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand. \url{https://www.otago.ac.nz/density/}.

}

\seealso{ 

\code{\link{capthist}}, 
\code{\link{occasionKey}}, 
\code{\link{sightingPlot}} 

}

\examples{
demotrap <- make.grid()
tempcapt <- sim.capthist(demotrap, 
    popn = list(D = 5, buffer = 50), 
    detectpar = list(g0 = 0.15, sigma = 30))
plot(tempcapt, border = 10, rad = 3, tracks = TRUE, 
    lab1cap = TRUE, laboffset = 2.5)

## type = n.per.cluster

## generate some captures
testregion <- data.frame(x = c(0,2000,2000,0),
    y = c(0,0,2000,2000))
popn <- sim.popn (D = 10, core = testregion, buffer = 0,
    model2D = "hills", details = list(hills = c(-2,3)))
t1 <- make.grid(nx = 1, ny = 1)
t1.100 <- make.systematic (cluster = t1, spacing = 100,
    region = testregion)
capt <- sim.capthist(t1.100, popn = popn, noccasions = 1)

## now plot captures ...
temp <- plot(capt, title = "Individuals per cluster",
    type = "n.per.cluster", hidetraps = FALSE,
    gridlines = FALSE, cappar = list(cex = 1.5))

if (interactive()) {
    ## add legend; click on map to place top left corner
    legend (locator(1), pch = 21, pt.bg = temp$colour,
        pt.cex = 1.3, legend = temp$legend, cex = 0.8)
}

\dontrun{

## try varying individual colours - requires RColorBrewer
library(RColorBrewer)
plot(infraCH[[2]], icolours = brewer.pal(12, "Set3"), tracks = TRUE,
    bg = "black", cappar = list(cex = 2), border = 10, rad = 2,
    gridlines = FALSE)

## generate telemetry data
te <- make.telemetry()
tr <- make.grid(detector = "proximity")
totalpop <- sim.popn(tr, D = 20, buffer = 100)
tepop <- subset(totalpop, runif(nrow(totalpop)) < 0.05)
teCH <- sim.capthist(te, popn = tepop, renumber=FALSE, detectfn = "HHN",
    detectpar = list(lambda0 = 3, sigma = 25))
plot(teCH, type = 'telemetry', tracks = TRUE)

## simple "centres" example
## polygon data require 'hazard' detection function 14:19
CH <- sim.capthist(make.poly(), nocc = 20, detectfn = 'HHN', 
    detectpar = list(lambda0 = 1, sigma = 10))
plot(CH, cappar = list(col = 'orange'), varycol = FALSE, border = 10)
plot(CH, type = 'centres', add = TRUE, rad = 0)
}

}

\keyword{ hplot }
