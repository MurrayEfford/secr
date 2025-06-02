\name{Rsurface}
\alias{Rsurface}

\title{ Smoothed Resource Surface}

\description{
  
  Creates a smoothed resource surface from a covariate of a
  mask. Smoothing entails summing the value in each pixel weighted by a
  detection kernel centred on the focal pixel. The detection kernel
  represents home-range utilization with spatial scale sigma. The
  resulting surface is equivalent to the denominator used by Royle et
  al. (2013) to normalize site-specific detection.
  
}

\usage{

Rsurface(mask, sigma, usecov = NULL, alpha2 = 1, detectfn = 'HHN', z = 1,
 inverse = FALSE, scale = TRUE)

}

\arguments{
  
  \item{mask}{secr habitat mask object (single-session)}
  \item{sigma}{numeric spatial scale of home range model}
  \item{alpha2}{numeric coefficient of spatial covariate}
  \item{usecov}{character name of resource covariate}
  \item{detectfn}{integer or character code for detection function}
  \item{z}{numeric shape parameter of home range model}
  \item{inverse}{logical; if TRUE the reciprocal of smoothed resource is returned}
  \item{scale}{logical; not used}

}

\details{

\code{detectfn} may be uniform (\sQuote{UN}) or one of the
cumulative-hazard functions (\sQuote{HHN}, \sQuote{HHR}, \sQuote{HEX},
\sQuote{HAN}, \sQuote{HCG}) (or integer codes 4, 14:18; see
\link{detectfn}).

The default \sQuote{HHN} corresponds to a halfnormal function on the hazard
scale, or a bivariate circular normal home range.

If \code{usecov} is not named then it takes the value 1.0 for all points
on the mask and zero otherwise.

The Rsurface can be used implicitly to normalize detection probability when
fitting a model with detector-specific covariate equal to
\code{usecov} (see \link{details}, but the process is intricate and not
fully documented).

}

\value{
  
An object with class c(`Rsurface', `mask', `data.frame') and covariate `Resource'
(other covariates are retained from the input mask). The attribute
`scale' is 1.0 if \code{scale = FALSE}; otherwise it is the average of the
resource over the masked area.

}

\references{

Royle, J. A., Chandler, R. B., Sun, C. C. and Fuller, A. K. (2013)
Integrating resource selection information with spatial
capture--recapture. \emph{Methods in Ecology and Evolution} \bold{4},
520--530.

}

\seealso{
  
  \code{\link{mask}}, \code{\link{plot.Rsurface}},
  \code{\link{spotHeight}}, \link{details}
  
}

\note{
  
Consider a focal pixel \bold{s} and another point in the habitat mask
\bold{x}, with distance \eqn{d = |x - s|}{d
= |x-s|}. Weights are given by a kernel \eqn{f(d)}. Typically the kernel
will be halfnormal \eqn{f(d) = \exp(-d^2/(2\sigma^2))}{f(d) =
exp(-d^2/(2 \sigma^2))} (detectfn = \sQuote{HHN}) or exponential
\eqn{f(d) = \exp(-d/\sigma)}{exp(-d/\sigma)} (detectfn =
\sQuote{HEX}) (see \link{detectfn} for other possibilities).

If \eqn{z(x)} represents the covariate value at point
\bold{x}, the summed resource availability at \bold{s} is given by
\deqn{R(s) = \sum_x f(d)\, \exp(\alpha_2 \,
  z(x)).}{R(s) = \sum f(d) exp(\alpha2 z(x)).}
This corresponds to the denominator of eqn 4 in Royle et al. (2013).

By default, the numerical values reported by \code{Rsurface} are not raw
\eqn{R} values. If \code{scale = TRUE}, values are standardized by
dividing by the mean: \eqn{R'(s) = R(s) /
  (\sum_s R(s) / n)}{R'(s) = R(s) / (\sum R(s) / n)} where \eqn{n} is the number of
pixels. Values of \eqn{R'(s)}{R'(s)} are centred on 1.0.

If \code{inverse = TRUE}, the numeric values are \eqn{1 /
  R'(s)}{1 / R'(s)}
or \eqn{1 / R(s)}{1 / R(s)} as determined by \code{scale}.

}

\examples{

\donttest{

## create binary covariate (0 outside habitat)
msk <- make.mask(traps(possumCH), buffer = 800)
covariates(msk) <- data.frame(z = as.numeric(pointsInPolygon
    (msk,possumarea)))

## derive and plot "resource availability"
Rs <- Rsurface(msk, sigma = 100, usecov = 'z')
plot(Rs, plottype = 'contour', col = topo.colors(10))
lines(possumarea)

if (interactive()) {
    spotHeight(Rs, dec = 2)
}

}

}
