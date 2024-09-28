\name{make.tri}
\alias{make.tri}
% \alias{make.hex}
\alias{clip.hex}

\title{ Build Detector Array on Triangular or Hexagonal Grid }

\description{
  Construct an array of detectors on a triangular grid and optionally
  select a hexagonal subset of detectors.
}

\usage{

make.tri (nx = 10, ny = 12, spacing = 20, detector = "multi",
    originxy = c(0,0))

clip.hex (traps, side = 20, centre = c(50, 60*cos(pi/6)),
    fuzz = 1e-3, ID = "num", ...)

}

\arguments{
  \item{nx}{ number of columns of detectors }
  \item{ny}{ number of rows of detectors }
  \item{spacing}{ distance between detectors (x and y directions)}
  \item{detector}{ character value for detector type - "single", "multi" etc. }
  \item{originxy}{ vector origin for x-y coordinates }

  \item{traps}{ traps object }
  \item{side}{ length of hexagon side }
  \item{centre}{ x-y coordinates of hexagon centre }
  \item{fuzz}{ floating point fuzz value }
  \item{ID}{ character string to control row names }
  \item{...}{ other parameters passed to subset.traps (not used) }
}

\details{

\code{make.tri} generates coordinates for \code{nx.ny} traps at
separations \code{spacing}. The bottom-left (southwest) corner is at
\code{originxy}. Identifiers are numeric. See \code{\link{make.grid}}
for further explanation.

\code{clip.hex} clips a grid of detectors, retaining only
those within a bounding hexagon. Detectors are re-labelled according to
\code{ID} as follows:

% \code{make.hex} generates coordinates for a hexagonal array obtained by
% thinning an initial triangular array of \code{nx.ny} traps at
% separations \code{spacing}. The bottom-left (southwest) corner is at
% \code{originxy}. Identifiers are numeric. 


  \tabular{ll}{
     ID \tab Effect \cr
     NULL \tab no change \cr
     num \tab numeric sequence \cr
     alpha \tab letter for`shell'; number within shell \cr
  }

}

\value{

An object of class \code{traps} comprising a data frame of x- and
y-coordinates, the detector type ("single", "multi", or "proximity" etc.),
and possibly other attributes.

}

\seealso{ \code{\link{make.grid}}, \code{\link{detector}}
}

\note{ Several methods are provided for manipulating detector arrays - see \code{\link{traps}}.

}

\examples{

tri.grid <- make.tri(spacing = 10)
plot(tri.grid, border = 5)

hex <- clip.hex(tri.grid, side = 30, ID = "alpha")
plot (hex, add = TRUE, detpar = list(pch = 16, cex = 1.4),
    label = TRUE, offset = 2.5 )

}


