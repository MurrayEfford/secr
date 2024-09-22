\name{Internal}

\alias{boundarytoSF}
\alias{Dfn2}
\alias{esaPlotsecr}

\title{Internal Functions}

\description{ 

Various functions called internally by \pkg{secr} and not exported.

}

\usage{

boundarytoSF (poly)

Dfn2(designD, beta = NULL, ...)

esaPlotsecr (object, max.buffer = NULL, max.mask = NULL,
    thin = 0.1, poly = NULL, poly.habitat = TRUE, session = 1, plt = TRUE, 
    type = "density", add = FALSE, overlay = TRUE, conditional = FALSE, ...)

}

\arguments{

  \item{poly}{data to define one or more polygons}

  \item{designD}{dataframe of density design data (output from \code{\link{D.designdata}})}
  \item{beta}{numeric vector of beta values (see Details for NULL)}

  \item{object}{\code{secr} object output from \code{secr.fit} }
  \item{max.buffer}{ maximum width of buffer in metres }
  \item{max.mask}{ \code{\link{mask}} object }
  \item{thin}{ proportion of mask points to retain in plot and output}
  \item{poly.habitat}{logical as in \code{\link{make.mask}} }
  \item{session}{vector of session indices (used if \code{object} spans multiple sessions)}
  \item{plt}{ logical to plot results }
  \item{type}{ character; see \code{\link{esaPlot}}}
  \item{add}{ logical to add line to an existing plot }
  \item{overlay}{ logical; if TRUE then automatically \code{add = TRUE}
  for plots after the first }
  \item{conditional}{logical; see \code{\link{esaPlot}}}
  \item{\dots}{ other arguments not used (boumndarytoSF, Dfn2) or graphical 
  arguments passed to plot() and lines() (esaPlotsecr)}

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

  \code{Dfn2} is supplied automatically as 'details' argument Dfn in 
  \code{\link{secr.fit}} when the switch Dlambda is set to TRUE for the 
  multi-session trend reparameterization of density. \code{Dfn2} 
  uses beta = NULL to return the required number of density coefficients 
  (beta parameters) in the model.

  \code{esaPlotsecr} provides a wrapper for \code{\link{esaPlot}}
  that is called internally from esaPlot when it is presented with an
  secr object. Arguments of \code{esaPlot} such as \code{detectfn} are
  inferred from the fitted model.

  The optional `poly' argument of \code{esaPlotsecr} is a matrix of two columns 
  interpreted as the x and y coordinates of a bounding polygon.

  If \code{max.mask} is not specified then a maximal mask of type
  `trapbuffer' is constructed using \code{max.buffer} and the spacing of
  the mask in \code{object}. In this case, if \code{max.buffer} is not
  specified then it is set either to the width of the existing plot
  (\code{add = TRUE}) or to 10 x sigma-hat from the fitted model in
  \code{object} (\code{add = FALSE}).

}

\value{

\code{boundarytoSF} -- Spatial object of \pkg{sf} class sfc, containing 
a geometry set of type POLYGON or MULTIPOLYGON. NULL input results in NULL output.

\code{Dfn2} -- Vector of density values on the link scale, suitable for the internal array (mask x groups x sessions).

\code{esaPlotsecr} -- see esaPlot

}

\seealso{

\code{\link{pointsInPolygon}}, 
\href{https://www.otago.ac.nz/density/pdfs/secr-spatialdata.pdf}{secr-spatialdata.pdf},

\code{\link{predictDlambda}},
\href{https://www.otago.ac.nz/density/pdfs/secr-trend.pdf}{secr-trend.pdf},


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