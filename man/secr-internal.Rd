\name{Internal}

\alias{boundarytoSF}

\title{Internal Functions}

\description{ 

Functions called internally by \pkg{secr} and not exported.

}

\usage{

boundarytoSF (poly)

}

\arguments{

  \item{poly}{data to define one or more polygons}

}

\details{

  The function \code{boundarytoSF} converts various 
  possible polygon input formats to a standard form (sfc).
  
  Possible inputs are:
  
  \tabular{lll}{
  Input \tab From \tab Note \cr
  2-column matrix or dataframe \tab base R\cr
  SpatialPolygons \tab \pkg{sp} \cr
  SpatialPolygonsDataFrame \tab \pkg{sp} \cr
  SpatVector \tab \pkg{terra} \cr
  sf \tab \pkg{sf} \tab geometry type POLYGON or MULTIPOLYGON \cr
  sfc \tab \pkg{sf} \tab geometry type POLYGON or MULTIPOLYGON \cr
  }
  
  Matrix input defines a single polygon.

  }

\value{

\code{boundarytoSF} -- Spatial object of \pkg{sf} class sfc, containing 
a geometry set of type POLYGON or MULTIPOLYGON. NULL input results in NULL output.

}

\seealso{

\code{\link{pointsInPolygon}}, \href{https://www.otago.ac.nz/density/pdfs/secr-spatialdata.pdf}{secr-spatialdata.pdf},

\code{\link{chat.nk}}

}

\references{

Hijmans, R. J. (2022) terra: Spatial Data Analysis. R package version 1.5-14. https://rspatial.org/terra/

Pebesma, E. (2018) Simple features for R: standardized support for spatial vector data. \emph{The R Journal} \bold{10(1)}, 439--446.
  https://doi.org/10.32614/RJ-2018-009

Pebesma, E.J. and Bivand, R. S. (2005) Classes and methods for spatial data in R. \emph{R News} \bold{5(2)}, 9--13.
  https://cran.r-project.org/doc/Rnews/Rnews_2005-2.pdf.

}

\examples{

\dontrun{

poly <- cbind(x = c(0,6,6,0,0), y = c(0,0,6,6,0))  
secr:::boundarytoSF(poly)

}

}