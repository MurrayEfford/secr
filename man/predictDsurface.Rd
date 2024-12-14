\name{predictDsurface}
\alias{predictDsurface}

\title{ Predict Density Surface }

\description{
Predict density at each point on a raster mask from a fitted secr model.
}

\usage{
predictDsurface(object, mask = NULL, se.D = FALSE, cl.D = FALSE, alpha =
0.05, parameter = c('D', 'noneuc'))
}

\arguments{
  \item{object}{fitted secr object}
  \item{mask}{secr mask object}
  \item{se.D}{ logical for whether to compute prediction SE}
  \item{cl.D}{ logical for whether to compute confidence limits}
  \item{alpha}{ alpha level for 100(1 -- alpha)\% confidence intervals }
  \item{parameter}{ character for real parameter to predict}
}

\details{

  Predictions use the linear model for density on the link scale in the
  fitted secr model `object', or the fitted user-defined function, if
  that was specified in \code{secr.fit}.
  
  If `mask' is NULL then predictions are for the mask component of
  `object'.

  SE and confidence limits are computed only if specifically
  requested. They are not available for user-defined density functions.

  Density is adjusted automatically for the number of clusters in
  `mashed' models (see \code{\link{mash}}).

  If \code{object} is a relative density model then \code{\link{derivedDsurface}} 
  first computes the intercept.
  
}

\value{

  Object of class `Dsurface' inheriting from `mask'. Predicted densities
  are added to the covariate dataframe (attribute `covariates') as
  column(s) with prefix `D.' If the model uses multiple groups, multiple
  columns will be distinguished by the group name (e.g., "D.F" and
  "D.M"). If groups are not defined the column is named "D.0".

  For multi-session models the value is a multi-session mask.
  
  The pointwise prediction SE is saved as a covariate column prefixed
  `SE.' (or multiple columns if multiple groups). Confidence limits are
  likewise saved with prefixes `lcl.'  and `ucl.'.
  
}

\seealso{
  \code{\link{plot.Dsurface}}, 
  \code{\link{secr.fit}}, 
  \code{\link{predict.secr}}, 
  \code{\link{derivedDsurface}}
}

\examples{

## use canned possum model
shorePossums <- predictDsurface(possum.model.Ds)
par(mar = c(1,1,1,6))
plot(shorePossums, plottype = "shaded", polycol = "blue", border = 100)
plot(traps(possumCH), detpar = list(col = "black"), add = TRUE)
par(mar = c(5,4,4,2) + 0.1)  ## reset to default
## extract and summarise
summary(covariates(shorePossums))

\dontrun{

## extrapolate to a new mask; add covariate needed by model; plot
regionmask <- make.mask(traps(possumCH), buffer = 1000, spacing = 10,
    poly = possumremovalarea)
dts <- distancetotrap(regionmask, possumarea)
covariates(regionmask) <- data.frame(d.to.shore = dts)
regionPossums <- predictDsurface(possum.model.Ds, regionmask,
    se.D = TRUE, cl.D = TRUE)
par(mfrow = c(1,2), mar = c(1,1,1,6))
plot(regionPossums, plottype = "shaded", mesh = NA, breaks = 20)
plot(regionPossums, plottype = "contour", add = TRUE)
plot(regionPossums, covariate = "SE", plottype = "shaded",
    mesh = NA, breaks = 20)
plot(regionPossums, covariate = "SE", plottype = "contour",
    add = TRUE)

## confidence surfaces
plot(regionPossums, covariate = "lcl", breaks = seq(0,3,0.2),
    plottype = "shaded")
plot(regionPossums, covariate = "lcl", plottype = "contour",
    add = TRUE, levels = seq(0,2.7,0.2))
title("lower 95\% surface")
plot(regionPossums, covariate = "ucl", breaks=seq(0,3,0.2),
    plottype = "shaded")
plot(regionPossums, covariate = "ucl", plottype = "contour",
    add = TRUE, levels = seq(0,2.7,0.2))
title("upper 95\% surface")

## annotate with CI
par(mfrow = c(1,1))
plot(regionPossums, plottype = "shaded", mesh = NA, breaks = 20)
plot(traps(possumCH), add = TRUE, detpar = list(col = "black"))

if (interactive()) {
    spotHeight(regionPossums, dec = 1, pre = c("lcl","ucl"), cex = 0.8)
}

## perspective plot
pm <- plot(regionPossums, plottype = "persp", box = FALSE, zlim =
    c(0,3), phi=30, d = 5, col = "green", shade = 0.75, border = NA)
lines(trans3d (possumremovalarea$x, possumremovalarea$y,
     rep(1,nrow(possumremovalarea)), pmat = pm))

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)  ## reset to default

## compare estimates of region N
## grid cell area is 0.01 ha
sum(covariates(regionPossums)[,"D.0"]) * 0.01
region.N(possum.model.Ds, regionmask)

}

}

