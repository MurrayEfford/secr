\name{userdist}
\alias{userdist}
\alias{noneuc}

\title{ Non-Euclidean Distances }

\description{

Non-Euclidean distances have a variety of uses, some obscure. You
probably do not need them unless you have data from linear habitats,
covered in the forthcoming package \pkg{secrlinear}. On the other hand,
they open up some intriguing possibilities for the advanced user. The
key is to provide an appropriate value for the component `userdist' of
the \code{details} argument of \code{\link{secr.fit}}.
  
\code{details$userdist} is either a function to compute distances
between detectors and mask points, or a pre-computed matrix of such
distances. Pre-computing assumes the matrix is static (i.e. fixed and
not dependent on any estimated coefficients). The functions
\code{\link{edist}} and \code{\link{nedist}} are useful for computing
static matrices of Euclidean or non-Euclidean distances (the latter is 
useful when there are barriers to movement).

If \code{details$userdist} is a function then it should take the form

userdist(xy1, xy2, mask)

}

\arguments{
  
\item{xy1}{2-column matrix of x-y coordinates of \eqn{k} detectors}
\item{xy2}{2-column matrix of x-y coordinates of \eqn{m} mask points}
\item{mask}{habitat mask defining a non-Euclidean habitat geometry}
}

\section{Details}{

The matrix returned by the function must have exactly \eqn{k} rows and
\eqn{m} columns. The function name may be almost anything you like.

The non-Euclidean habitat geometry may or may not require access to
local density (D), local (mask) covariates, and the estimation of
additional coefficients (beta variables). In order that secr.fit can
assemble these data, there is a mechanism for the user to indicate
which, if any, variables are required: when called with no arguments the
function should return a character vector of variable names. These may
include covariates of `mask', the dynamically computed density `D`, and
a special real parameter `noneuc' for which one or more coefficients
will be fitted.

`noneuc' is like 'D' in that it may be modelled as a function of any mask
covariates, session, Session, x, y, etc. The actual meaning attributed
to `noneuc' depends entirely on how it is used inside the function.

The function may require no variables and not require estimation of
additional coefficients. This is the case for a simple linear geometry
as described in documentation for the package `secrlinear'.

\tabular{ll}{
  Value  \tab Interpretation \cr
  ''       \tab no covariates etc. required \cr
  'D'      \tab density at each mask point \cr
  'noneuc' \tab a multi-purpose real parameter \cr
  \tab defined for each mask point \cr
  c('D', 'noneuc') \tab both of the preceding \cr
  c('noneuc','habclass') \tab both noneuc and the mask covariate
  'habclass' \cr  
}

The last case does not estimate a coefficient for habclass, it merely
makes the raw value available to whatever algorithm you implement. 

The `xy2' and `mask' parameters of the userdist function overlap in
practice: xy1 and xy2 only define the points between which distances are
required, whereas mask is a carrier for any and all additional
information needed by the algorithm.

Full documentation of the \pkg{secr} capability for non-Euclidean
distances is in the separate document
\href{https://www.otago.ac.nz/density/pdfs/secr-noneuclidean.pdf}{secr-noneuclidean.pdf},
which includes example code for the analysis of Sutherland et
al. (2015).

}

\section{Compatibility}{
  
User-specified distances are compatible with some but not all features
of \pkg{secr}. Functions with a `userdist' argument are certainly
compatible, and others may be.

With a static userdist, \code{region.N} will generally not calculate  population size for a region other than the original mask. If you want to supply a new mask in the `region' argument, replace x$details$userdist with a distance matrix appropriate to the new mask, where `x' is the name of the fitted model.

User-specified distances cannot be used with polygon or transect
detectors.

When using \code{\link{sim.capthist}} to simulate detections of a new
population from \code{\link{sim.popn}} you must provide \code{userdist}
as a function rather than a matrix. This is because new animals are not
restricted to locations on the `mask' grid. 
  
}

\seealso{
  \code{\link{details}}, \code{\link{secr.fit}}, \code{\link{nedist}}
}

\references{

Sutherland, C., Fuller, A. K. and Royle, J. A. (2015) Modelling
non-Euclidean movement and landscape connectivity in highly structured
ecological networks. \emph{Methods in Ecology and Evolution} \bold{6}, 169--177.

}

\examples{

## see secr-noneuclidean.pdf

}

\keyword{ models }
