\name{mask}
\alias{mask}
\alias{mask-class}

\title{ Mask Object }

\description{
Encapsulate a habitat mask for spatially explicit capture--recapture. See also \href{https://www.otago.ac.nz/density/pdfs/secr-habitatmasks.pdf}{secr-habitatmasks.pdf}.

}

\details{

A habitat mask serves four main purposes in spatially explicit
capture--recapture. Firstly, it defines an outer limit to the area of
integration; habitat beyond the mask may be occupied, but animals there
should have negligible chance of being detected (see \code{\link{pdot}}
and below). Secondly, it distinguishes sites in the vicinity of the
detector array that are `habitat' (i.e. have the potential to be
occupied) from `non-habitat'. Thirdly, it discretizes continuous habitat
as a list of points. Each point is notionally associated with a cell
(pixel) of uniform density. Discretization allows the SECR likelihood to
be evaluated by summing over grid cells. Fourthly, the x-y coordinates
of the mask and any habitat covariates may be used to build spatial
models of density. For example, a continuous or categorical habitat
covariate `cover' measured at each point on the mask might be used in a
formula for density such as D \eqn{\sim}{~}cover.

In relation to the first purpose, the definition of `negligible' is
fluid. Any probability less than 0.001 seems OK in the sense of not
causing noticeable bias in density estimates, but this depends on the
shape of the detection function (fat-tailed functions such as `hazard
rate' are problematic). New tools for evaluating masks appeared in
\pkg{secr} 1.5 (\code{\link{mask.check}}, \code{\link{esaPlot}}), and
\code{\link{suggest.buffer}} automates selection of a buffer width.

Mask points are stored in a data frame with columns `x' and `y'. The
number of rows equals the number of points.


Possible mask attributes
\tabular{ll}{
Attribute \tab Description \cr
type \tab `traprect', `trapbuffer', `pdot', `polygon', `clusterrect', `clusterbuffer' (see make.mask) or `user' \cr
polygon \tab vertices of polygon defining habitat boundary, for type = `polygon' \cr
pdotmin \tab threshold of p.(X) for type = `pdot' \cr
covariates \tab dataframe of site-specific covariates \cr
meanSD \tab data frame with centroid (mean and SD) of x and y coordinates \cr
area \tab area (ha) of the grid cell associated with each point \cr
spacing \tab nominal spacing (metres) between adjacent points \cr
boundingbox \tab data frame of 4 rows, the vertices of the bounding box of all grid cells in the mask \cr
}

Attributes other than \code{covariates} are generated automatically by
\code{make.mask}. Type `user' refers to masks input from a text file
with \code{read.mask}.

A virtual S4 class `mask' is defined to allow the definition of a method
for the generic function \code{raster} from the \pkg{raster} package.

}

\note{ A habitat mask is needed by \code{secr.fit}, but one will be
generated automatically if none is provided. You should be aware of this
and check that the default settings (e.g. \code{buffer}) are
appropriate.  }

\seealso{ 
  \code{\link{make.mask}}, \code{\link{read.mask}},
  \code{\link{mask.check}}, \code{\link{esaPlot}},
  \code{\link{suggest.buffer}}, \code{\link{secr.fit}}
}

\keyword{ classes }

