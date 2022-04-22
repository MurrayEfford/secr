\name{Dsurface}
\alias{Dsurface}
\alias{Dsurface-class}
\alias{print.Dsurface}
\alias{summary.Dsurface}

\title{ Density Surfaces }

\description{

  S3 class for rasterized fitted density surfaces. A Dsurface is a type
  of `mask' with covariate(s) for the predicted density at each point.

}

\usage{

\method{print}{Dsurface}(x, scale = 1, ...)
\method{summary}{Dsurface}(object, scale = 1, ...)

}

\arguments{

  \item{x, object}{Dsurface object to display}
  \item{scale}{numeric multiplier for density}
  \item{...}{other arguments passed to print method for data frames or
  summary method for masks}

}

\details{
    
A Dsurface will usually have been constructed with \code{\link{predictDsurface}}.

The `scale' argument may be used to change the units of density from the default
(animals / hectare) to animals / km^2 (scale = 100) or animals / 100km^2
(scale = 10000).

A virtual S4 class `Dsurface' is defined to allow the definition of a method for the generic function \code{raster} from the \pkg{raster} package.

}

\seealso{
  
\code{\link{predictDsurface}}, \code{\link{plot.Dsurface}}

}

\keyword{ classes }
