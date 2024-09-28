\name{fxTotal}
\alias{fxTotal}
\alias{fxTotal.secr}

\title{Activity Centres of Detected and Undetected Animals}

\description{

  The summed probability densities of both observed and
  unobserved individuals are computed for a fitted model and dataset.

  Function \code{fx.total} was replaced by method \code{fxTotal} in \pkg{secr} 5.0.0.
}

\usage{

\method{fxTotal}{secr}(object, sessnum = 1, mask = NULL, ncores = NULL, ...)

}

\arguments{
  \item{object}{ a fitted secr model }
  
  \item{sessnum}{ session number if \code{object$capthist} spans
    multiple sessions }
  
  \item{mask}{ x- and y- coordinates of points at which density will be
    computed}
    
  \item{ncores}{ integer number of threads to be used for parallel processing}
    
  \item{\dots}{ other arguments passed to \code{detectpar} and thence
    to \code{predict.secr}}
    
  }
  
\details{

This function calls \code{\link{fxi}} for each detected animal and
overlays the results to obtain a summed probability density surface D.fx
for the locations of the home-range centres of detected individuals.

A separate calculation using \code{\link{pdot}} provides the expected
spatial distribution of undetected animals, as another density
surface: crudely, D.nc(X) = D(X) * ( 1 -- pdot(X)).

The pointwise sum of the two surfaces is sometimes used to represent the
spatial distrbution of the population, but see Notes.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

}

\value{
An object of class `Dsurface' (a variety of mask) with a `covariates'
attribute that is a dataframe with columns --
\item{D.fx }{sum of \code{\link{fxi}} over all detected individuals}
\item{D.nc }{expected density of undetected (`not caught') individuals}
\item{D.sum }{sum of D.fx and D.nc}

All densities are in animals per hectare (the `scale' argument of
\code{\link{plot.Dsurface}} allows the units to be varied later).

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

}

\note{

The surface D.sum represents what is known from the data about a
specific realisation of the spatial point process for home range
centres: varying the intensity of sampling will change its shape. It is
not an unbiased estimate of a biologically meaningful density
surface. The surface will always tend to lack relief towards the edge of
a habitat mask where the main or only contribution is from D.nc.

}

\seealso{
  
  \code{\link{fxi}}, \code{\link{fxiContour}}, \code{\link{pdot}}
  
}

\examples{

\dontrun{

tmp <- fxTotal(secrdemo.0)

## to plot we must name one of the covariates:
## the Dsurface default 'D.0' causes an error 

plot(tmp, covariate = 'D.sum', col = terrain.colors(16),
   plottype = 'shaded')
plot(tmp, covariate = 'D.sum', col = 'white', add = TRUE,
   plottype = 'contour')
if (interactive()) {
    spotHeight(tmp, prefix = 'D.sum')
}

fxsurface <- fxTotal(ovenbird.model.D, sessnum = 3)
plot(fxsurface, covariate = 'D.sum')

}


}


