\name{summary.mask}
\alias{summary.mask}
\alias{print.summary.mask}

\title{ Summarise Habitat Mask }

\description{
  Concise summary of a \code{mask} object.
}
\usage{
\method{summary}{mask}(object, ...)
\method{print}{summary.mask}(x, ...)
}

\arguments{
  \item{object}{ \code{mask} object }
  \item{x}{ \code{summary.mask} object }
  \item{\dots}{ other arguments (not used) }
}

\details{
The bounding box is the smallest rectangular area with edges parallel to the x- and y-axes that contains all points and their associated grid cells. A print method is provided for objects of class \code{summary.mask}.
}

\value{

  Object of class `summary.mask', a list with components
  \item{detector }{character string for detector type ("single","multi","proximity")}
  \item{type }{mask type ("traprect", "trapbuffer", "pdot", "polygon", "user", "subset")}
  \item{nmaskpoints }{number of points in mask}
  \item{xrange }{range of x coordinates}
  \item{yrange }{range of y coordinates}
  \item{meanSD }{dataframe with mean and SD of x, y, and each covariate }
  \item{spacing }{nominal spacing of points}
  \item{cellarea }{area (ha) of grid cell associated with each point}
  \item{bounding box}{dataframe with x-y coordinates for vertices of bounding box}
  \item{covar }{summary of each covariate}
}

\seealso{ \code{\link{mask}} }

\examples{
tempmask <- make.mask(make.grid())
## left to right gradient
covariates (tempmask) <- data.frame(x = tempmask$x)
summary(tempmask)
}

\keyword{ models }

