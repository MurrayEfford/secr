\name{esa.plot.secr}
\alias{esa.plot.secr}

\title{

Mask Buffer Diagnostic Plot (internal)

}

\description{
  
Internal function used to plot effective sampling area (Borchers and
Efford 2008) as a function of increasing buffer width given an `secr'
object

}

\usage{

esa.plot.secr (object, max.buffer = NULL, max.mask = NULL,
    thin = 0.1, poly = NULL, poly.habitat = TRUE, session = 1, plt = TRUE, 
    type = "density", add = FALSE, overlay = TRUE, conditional = FALSE, ...)

}

\arguments{
  \item{object}{\code{secr} object output from \code{secr.fit} }
  \item{max.buffer}{ maximum width of buffer in metres }
  \item{max.mask}{ \code{\link{mask}} object }
  \item{thin}{ proportion of mask points to retain in plot and output}
  \item{poly}{ matrix of two columns interpreted as the x and y
    coordinates of a bounding polygon (optional) }
  \item{poly.habitat}{logical as in \code{\link{make.mask}} }
  \item{session}{vector of session indices (used if \code{object} spans multiple sessions)}
  \item{plt}{ logical to plot results }
  \item{type}{ character; see \code{\link{esa.plot}}}
  \item{add}{ logical to add line to an existing plot }
  \item{overlay}{ logical; if TRUE then automatically \code{add = TRUE}
  for plots after the first }
  \item{conditional}{logical; see \code{\link{esa.plot}}}
  \item{\dots}{ graphical arguments passed to plot() and lines()}
}

\details{

  \code{esa.plot.secr} provides a wrapper for \code{\link{esa.plot}}
  that is called internally from esa.plot when it is presented with an
  secr object. Arguments of \code{esa.plot} such as \code{detectfn} are
  inferred from the fitted model.

  If \code{max.mask} is not specified then a maximal mask of type
  `trapbuffer' is constructed using \code{max.buffer} and the spacing of
  the mask in \code{object}. In this case, if \code{max.buffer} is not
  specified then it is set either to the width of the existing plot
  (\code{add = TRUE}) or to 10 x sigma-hat from the fitted model in
  \code{object} (\code{add = FALSE}).
  
}

\value{
  see \code{esa.plot}
}

\seealso{ 

\code{\link{esa.plot}}, \code{\link{mask}},
\code{\link{pdot}}, \code{\link{make.mask}}, \code{\link{mask.check}},
\code{\link{Detection functions}} 

}

\keyword{ hplot }
