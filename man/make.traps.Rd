\name{make.traps}
\alias{make.grid}
\alias{make.circle}
\alias{make.poly}
\alias{make.transect}
\alias{make.telemetry}

\title{ Build Detector Array }

\description{
  Construct a rectangular array of detectors (trapping grid) or a circle of detectors or a polygonal search area.
}

\usage{
make.grid(nx = 6, ny = 6, spacex = 20, spacey = spacex, spacing = NULL,
    detector = "multi", originxy = c(0,0), hollow = FALSE, ID = "alphay", 
    leadingzero = TRUE, markocc = NULL)

make.circle (n = 20, radius = 100, spacing = NULL, 
    detector = "multi", originxy = c(0,0), IDclockwise = TRUE, 
    leadingzero = TRUE, markocc = NULL)

make.poly (polylist = NULL, x = c(-50,-50,50,50),
    y = c(-50,50,50,-50), exclusive = FALSE, verify = TRUE)

make.transect (transectlist = NULL, x = c(-50,-50,50,50),
    y = c(-50,50,50,-50), exclusive = FALSE)

make.telemetry (xy = c(0,0)) 

}

\arguments{
  \item{nx}{ number of columns of detectors }
  \item{ny}{ number of rows of detectors }
  \item{spacex}{ distance between detectors in `x' direction (nominally in metres)}
  \item{spacey}{ distance between detectors in `y' direction (nominally in metres)}
  \item{spacing}{ distance between detectors (x and y directions)}
  \item{detector}{ character value for detector type - "single", "multi" etc. }
  \item{originxy}{ vector origin for x-y coordinates }
  \item{hollow}{ logical for hollow grid }
  \item{ID}{ character string to control row names }
  \item{leadingzero}{logical; if TRUE numeric rownames are padded with leading zeros}
  \item{markocc}{integer vector of marking or sighting codes; see \code{\link{markocc}}}
  \item{n}{ number of detectors }
  \item{radius}{ radius of circle  (nominally in metres) }
  \item{IDclockwise}{ logical for numbering of detectors }
  \item{polylist}{ list of dataframes with coordinates for polygons }
  \item{transectlist}{ list of dataframes with coordinates for transects }
  \item{x}{ x coordinates of vertices }
  \item{y}{ y coordinates of vertices }
  \item{exclusive}{ logical; if TRUE animal can be detected only once
  per occasion }
  \item{verify}{ logical if TRUE then the resulting traps object is
    checked with \code{\link{verify}}}
  \item{xy}{vector with coordinates of arbitrary point (e.g., centroid of fixes)}
    
}

\details{

\code{make.grid} generates coordinates for \code{nx.ny} traps at
separations \code{spacex} and \code{spacey}. If \code{spacing} is
specified it replaces both \code{spacex} and \code{spacey}. The
bottom-left (southwest) corner is at \code{originxy}. For a hollow grid,
only detectors on the perimeter are retained. By default, identifiers
are constructed from a letter code for grid rows and an integer value
for grid columns ("A1", "A2",...). `Hollow' grids are always numbered
clockwise in sequence from the bottom-left corner. Other values of
\code{ID} have the following effects:

  \tabular{ll}{
     ID \tab  Effect \cr
     numx \tab  column-dominant numeric sequence \cr
     numy \tab  row-dominant numeric sequence \cr
     numxb \tab  column-dominant boustrophedonical numeric sequence (try it!) \cr
     numyb \tab  row-dominant boustrophedonical numeric sequence \cr
     alphax \tab  column-dominant alphanumeric \cr
     alphay \tab  row-dominant alphanumeric \cr
     xy \tab  combine column (x) and row(y) numbers \cr     
  }

  `xy' adds leading zeros as needed to give a string of constant length with
  no blanks.
  
\code{make.circle} generates coordinates for n traps in a circle centred
on \code{originxy}. If \code{spacing} is specified then it overrides the
\code{radius} setting; the radius is adjusted to provide the requested
straightline distance between adjacent detectors. Traps are numbered
from the trap due east of the origin, either clockwise or anticlockwise
as set by \code{IDclockwise}. 


Polygon vertices may be specified with \code{x} and \code{y} in the case
of a single polygon, or as \code{polylist} for one or more polygons. Each
component of \code{polylist} is a dataframe with columns `x' and `y'.
\code{polylist} takes precedence. \code{make.poly} automatically closes
the polygon by repeating the first vertex if the first and last vertices
differ.

Transects are defined by a sequence of vertices as for polygons, except 
that they are not closed.

\code{make.telemetry} builds a simple traps object for the 'telemetry' detector type. The attribute 'telemetrytype' is set to "independent". 

Specialised functions for arrays using a triangular grid are described
separately (\code{\link{make.tri}},  \code{\link{clip.hex}}).

}

\value{

An object of class \code{traps} comprising a data frame of x- and
y-coordinates, the detector type ("single", "multi", or "proximity" etc.),
and possibly other attributes.

}

\references{ 

Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially explicit
capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand. \url{https://www.otago.ac.nz/density/}.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

}

\seealso{
  \code{\link{read.traps}},\code{\link{detector}},
  \code{\link{trap.builder}},\code{\link{make.systematic}},  
  \code{\link{print.traps}}, \code{\link{plot.traps}},
  \code{\link{traps}}, \code{\link{make.tri}},
  \code{\link{addTelemetry}}
}

\note{ Several methods are provided for manipulating detector arrays - see \code{\link{traps}}.

}

\examples{
demo.traps <- make.grid()
plot(demo.traps)

## compare numbering schemes
par (mfrow = c(2,4), mar = c(1,1,1,1), xpd = TRUE)
for (id in c("numx", "numy", "alphax", "alphay", "numxb", 
    "numyb"))
{
    temptrap <- make.grid(nx = 7, ny = 5, ID = id)
    plot (temptrap, border = 10, label = TRUE, offset = 7, 
        gridl = FALSE)
}

temptrap <- make.grid(nx = 7, ny = 5, hollow = TRUE)
plot (temptrap, border = 10, label = TRUE, gridl = FALSE)

plot(make.circle(n = 20, spacing = 30), label = TRUE, offset = 9)
summary(make.circle(n = 20, spacing = 30))

## jitter locations randomly within grid square
## and plot over `mask'
# see also ?gridCells

tr0 <- tr <- make.grid(nx = 7, ny = 7, spacing = 30)
tr[] <- jitter(unlist(tr), amount = spacing(tr)/2)
plot(as.mask(tr0), dots = FALSE, mesh = 'white')
plot(tr, add = TRUE)

}

\keyword{ datagen }

