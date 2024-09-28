\name{smooths}
\alias{smooths}

\title{ Smooth Terms in SECR Models }

\description{

  From version 2.9.0, the model formulae provided to \code{secr.fit} may
  include smooth terms as specified for the \pkg{mgcv} function `gam',
  with some restrictions. Smooth terms may be used for both density and
  detection parameters.

}

\section{Background}{

  Regression splines are a very flexible way to represent non-linear
  responses in generalized additive models (e.g., \pkg{mgcv}, Wood
  2006). Borchers and Kidney (2014) have shown how they may be used
  to model 2-dimensional trend in density in \pkg{secrgam}, an \R
  package that extends \pkg{secr}. Their approach is to use \pkg{mgcv}
  to construct regression spline basis functions from mask x- and
  y-coordinates, and possibly additional mask covariates, and to pass
  these as covariates to \pkg{secr}. The idea of using \pkg{mgcv} to
  construct the basis functions is applied within \pkg{secr} from version 2.9.

  Smooth semi-parametric responses are also useful for modelling
  variation in detection parameters such as g0 and sigma over time, or
  in response to individual or detector-level covariates, when (1) a
  linear or other parametric response is arbitrary or implausible, and
  (2) sampling spans a range of times or levels of the covariate(s).

  For a concrete example, consider a population sampled monthly for a
  year (i.e., 12 `sessions'). If home range size varies seasonally then
  the parameter sigma may vary in a more-or-less sinusoidal fashion. A
  linear trend is obviously inadequate, and a quadratic is not much
  better. However, a sine curve is hard to fit (we would need to
  estimate its phase, amplitude, mean and spatial scale) and assumes the
  increase and decrease phases are equally steep. An extreme solution is
  to treat month as a factor and estimate a separate parameter for each
  level (month). A smooth (semi-parametric) curve may capture the main
  features of seasonal variation with fewer parameters.
}

\details{

  The specification of smooth terms is explained in
  \code{\link[mgcv]{formula.gam}}. Only a subset of options are relevant
  to `secr'. Penalized splines are not available. The smooth function
  may be `s' or `te'.
  
  The `wiggliness' of the curve is controlled by the argument k,
  which in this implementation is set by the user. The argument `fx'
  should be set to TRUE.

%  Polynomial regression is also allowed. The function \code{poly}
%  specifies a polynomial regression using orthogonal polynomials (e.g.,
%  \code{g0 ~ poly(T, 2) for a quadratic}). This is both convenient and
%  efficient (unlike the usual terms in x, x^2 etc., orthogonal
%  polynomials are uncorrelated). Note that \code{poly} may include more
%  than one predictor, but in this case the degree argument must be
%  named (e.g., \code{poly(x, y, degree = 2)}).

See also the example in \href{https://www.otago.ac.nz/density/pdfs/secr-densitysurfaces.pdf}{secr-densitysurfaces.pdf}.

}

\note{

  Setting the intercept of a smooth to zero is not a canned option in
  \pkg{mgcv}, and is not offered in \pkg{secr}. It may be achieved by
  placing a knot at zero and hacking the matrix of basis functions to
  drop the corresponding column, plus some more jiggling.
  %\href{https://r.789695.n4.nabble.com/Use-pcls-in-quot-mgcv-quot-package-to-achieve-constrained-cubic-spline-td4660966.html}{jiggling}. NOT FOUND. TO BE CHECKED

}

\section{Warning}{

  This implementation of smooth models results in large fitted objects, on
  account of the need to store setup information from \pkg{mgcv}. It is
  also vulnerable to future changes in \pkg{mgcv}.

  Expect that the implementation will change in later versions of
  \pkg{secr}, and that smooth models fitted in the this version will not
  necessarily be compatible with \code{predict} and
  \code{predictDsurface} in later versions.
  
}

\examples{

\dontrun{

## smooth density surface
possum.model.sxy <- secr.fit(possumCH, mask = possummask,
    model = D ~ s(x,y, k = 6, fx = TRUE), trace = FALSE)
fittedsurface <- predictDsurface(possum.model.sxy)
par(mar = c(1,1,1,6))
plot(fittedsurface)
plot(fittedsurface, plottype = 'contour', add = TRUE)
par(mar = c(5,4,4,2) + 0.1)  ## reset to default

## Now try smooth on g0

## For the smooth we use 'Session' which is coded numerically (0:4)
## rather than the factor 'session' ('2005', '2006', '2007', '2008',
## '2009')

ovenbird.model.g0 <- secr.fit(ovenCH, mask = ovenmask,
    model = g0 ~ session, trace = FALSE)
ovenbird.model.sg0 <- secr.fit(ovenCH, mask = ovenmask,
    model = g0 ~ s(Session, k = 3, fx = TRUE), trace = FALSE)

AIC(ovenbird.model.g0, ovenbird.model.sg0)

## Or over occasions within a session...

fit.sT3 <- secr.fit(captdata, model = g0 ~ s(T, k = 3, fx = TRUE),
    trace = FALSE)
pred <- predict(fit.sT3, newdata = data.frame(T = 0:4))

plot(sapply(pred, '[', 'g0', 'estimate'))

}

}

\references{ 

Borchers, D. L. and Kidney, D. (2014) \emph{Flexible density surface estimation for spatially explicit capture--recapture surveys}. Technical Report, University of St Andrews.

Wood, S. N. (2006) \emph{Generalized Additive Models: An Introduction with
R.} Chapman and Hall/CRC.

}

\seealso{
  
  \code{\link[mgcv]{formula.gam}} % , \code{\link[stats]{poly}}
  
}


