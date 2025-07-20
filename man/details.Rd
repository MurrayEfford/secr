\name{details}
\alias{details}
\alias{fixedbeta}
\alias{LLonly}
\alias{param}
\alias{miscparm}
\alias{maxdistance}
\alias{saveprogress}

\title{ Detail Specification for secr.fit }

\description{

The function \code{secr.fit} allows many options. Some of these are used
infrequently and have been bundled as a single argument \code{details}
to simplify the documentation. They are described here. Some components (param, 
relativeD) are primarily used to record the type of model fitted that is an 
indirect result of other settings.

}

\section{Detail components}{

\code{details$autoini} specifies the session number from which to compute starting 
values (multi-session data only; default 1). From 4.0.0, the character value `all' 
first forms a single-session capthist using \code{\link{join}}(); this may be slow or not 
work at all (especially with telemetry data).

\code{details$centred} = TRUE causes coordinates of both traps and mask
to be centred on the centroid of the traps, computed separately for each
session in the case of multi-session data. This may be necessary to
overcome numerical problems when x- or y-coordinates are large
numbers. The default is not to centre coordinates.

\code{details$chat} optionally specifies the overdispersion 
of unmarked sightings Tu and unidentified marked sightings Tm. It is used only 
for mark-resight models, and is usually computed within \code{secr.fit} 
(\code{details$nsim > 0}), but may be provided by the user. For a single session
`chat` is a vector of length 2; for multiple sessions it is a 2-column matrix. 

\code{details$chatonly} = TRUE used with \code{details$nsim > 0} causes the 
overdispersion statistics for sighting counts Tu and Tm to be estimated and 
returned as a vector or 2-column matrix (multi-session models), with no further
model fitting.

\code{details$contrasts} may be used to specify the coding of factor predictors. 
The value should be suitable for the 'contrasts.arg' argument of 
\code{\link[stats]{model.matrix}}. See `Trend across sessions' in
\href{https://www.otago.ac.nz/density/pdfs/secr-multisession.pdf}{secr-multisession.pdf} 
for an example.

\code{details$convexpolygon} may be set to FALSE for searches of non-convex polygons. 
This is slower than the default which requires poygons to be convex east-west 
(\href{https://www.otago.ac.nz/density/pdfs/secr-polygondetectors.pdf}{secr-polygondetectors.pdf}).

\code{details$debug} is an integer code used to control the printing of intermediate 
values (1,2) and to switch on the R code browser (3). In ordinary use it should not be 
changed from the default (0).

\code{details$Dfn} is a function for reparameterizing density models; this is set 
internally when Dlambda = TRUE. Exotic variations may be specified directly by 
the user when Dlambda = FALSE. The defaults (Dfn = NULL, Dlambda = FALSE) leave
the original density model unchanged. Note there is no connection to userDfn 
(except that the two are incompatible).

\code{Dlambda} if TRUE causes reparameterization of density as the 
session-on-session finite rate of increase \eqn{lambda}. Details at 
(\href{https://www.otago.ac.nz/density/pdfs/secr-trend.pdf}{secr-trend.pdf}).

\code{details$distribution} specifies the distribution of the number of
individuals detected \eqn{n}; this may be conditional on the number in the
masked area ("binomial") or unconditional ("poisson").
\code{distribution} affects the sampling variance of the estimated
density. The default is "poisson". The component `distribution' may also
take a numeric value larger than nrow(capthist), rather than "binomial"
or "poisson". The likelihood then treats n as a binomial draw from a
superpopulation of this size, with consequences for the variance of
density estimates. This can help to reconcile MLE with Bayesian
estimates using data augmentation.

\code{details$externalpdot} names a mask covariate that is used as an offset for 
relative density. This can be useful in a two-phase study when animals are tagged 
in phase one and sampled in phase two, with no further tagging (Efford in prep.). 
The covariate may differ from \eqn{p.(x)} by a constant factor. 

\code{details$fastproximity} controls special handling of data from binary 
proximity and count detectors. If TRUE and other conditions are met (no temporal 
variation or groups) then a multi-occasion capthist is automatically reduced to 
a count for a single occasion and further compressed by storing only non-zero 
counts, which can greatly speed up computation of the likelihood (default TRUE).

\code{details$fixedbeta} may be used to fix values of beta
parameters. It should be a numeric vector of length equal to the total
number of beta parameters (coefficients) in the model. Parameters to be
estimated are indicated by NA. Other elements should be valid values on
the link scale and will be substituted during likelihood
maximisation. Check the order of beta parameters in a previously fitted
model.

\code{details$grain} sets the grain argument for multithreading in RcppParallel 
parallelFor (default 1). \code{details$grain = 0} suppresses multithreading 
(equivalent to \code{ncores = 1}).

\code{details$hessian} is a character string controlling the computation
of the Hessian matrix from which variances and covariances are obtained.
Options are "none" (no variances), "auto" (the default) or "fdhess" (use
the function fdHess in \pkg{nlme}).  If "auto" then the Hessian from the
optimisation function is used.

\code{details$ignoreusage} = TRUE causes the function to ignore
usage (varying effort) information in the traps component. The default
(\code{details$ignoreusage} = FALSE) is to include usage in the model.

\code{details$intwidth2} controls the half-width of the interval
searched by optimise() for the maximum likelihood when there is a single
parameter. Default 0.8 sets the search interval to \eqn{(0.2s, 1.8s)} where \eqn{s}
is the `start' value.

\code{details$knownmarks} = FALSE causes secr.fit to fit a zero-truncated 
sightings-only model that implicitly estimates the number of marked individuals,
rather than inferring it from the number of rows in the capthist object.

\code{details$LLonly} = TRUE causes the function to returns a single
evaluation of the log likelihood at the `start' values.

\code{details$maxdistance} sets a limit to the centroid-to-mask distances 
considered. The centroid is the geometric mean of detection locations for each 
individual. If no limit is specified then summation is over all mask points. 
Specifying \code{maxdistance} can speed up computation; it is up to the user to 
select a limit that is large enough not to affect the likelihood (\eqn{5\sigma}?).

\code{details$miscparm} (default NULL) is an optional numeric vector of
starting values for additional parameters used in a user-supplied
distance function (see `userdist' below). If the vector has a names
attribute then the names will be used for the corresponding coefficients
(`beta' parameters) which will otherwise be named `miscparm1',
miscparm2' etc. These parameters are constant across each model and do
not appear in the model formula, but are estimated along with other
coefficients when the likelihood is maximised. Any transformation (link
function) etc. is handled by the user in the userdist function. The
coefficients appear in the output from \code{\link{coef.secr}} and
\code{\link{vcov.secr}}, but not \code{predict.secr}.

\code{details$newdetector} specifies a detector type to use for this fit, 
replacing the previous \code{detector(traps(capthist))}. The value may be 
a vector (one value per occasion) or for multi-session data, a list of vectors. 
A scalar value (e.g. "proximity") is otherwise used for all occasions and sessions. 
The true detector type is usually known and will be specified in the 'traps' 
attribute; \code{newdetector} is useful in simulation studies that examine the 
effect of misspecification. The capthist component of the output from secr.fit 
has the new type.

\code{details$nsim} specifies the number of replicate simulations to 
perform to estimate the overdispersion statistics for the sighting counts 
Tu and Tm. See also \code{details$chat} and \code{details$chatonly}.

\code{details$param} chooses between various parameterisations of the
SECR model. The default \code{details$param = 0} is the formulation in
Borchers and Efford (2008) and later papers.

\code{details$param = 1} was once used to select the Gardner & Royle parameterisation of
the detection model (\emph{p0}, \eqn{\sigma}; Gardner et al. 2009) when
the detector type is `multi'. This parameterisation was discontinued in 2.10.0.

\code{details$param = 2} selects parameterisation in terms of
(\eqn{esa(g_0, \sigma)}, \eqn{\sigma}) (Efford and Mowat 2014).

\code{details$param = 3} selects parameterisation in terms of
(\eqn{a_0(\lambda_0, \sigma)}, \eqn{\sigma}) (Efford and Mowat 2014). This
parameterization is used automatically if a0 appears in the model (e.g.,
a0 ~ 1).

\code{details$param = 4} selects parameterisation of sigma in terms of
the coefficient sigmak and constant c (sigma = sigmak /
  D^0.5 + c) (Efford et al. 2016). If c is not included explicitly in
  the model (e.g., c ~ 1) then it is set to zero. This
parameterization is used automatically if sigmak appears in the model (e.g.,
sigmak ~ 1)

\code{details$param = 5} combines parameterisations (3) and (4) (first
compute sigma from D, then compute lambda0 from sigma).

\code{details$relativeD} indicates that a density model has been fitted conditional 
on \eqn{n} and hence describing relative density instead of absolute density. 
The intercept of the submodel for D as a function of covariates is 
fixed at zero for a log link or 1 for the identity link. 
Not set by the user in secr >= 5.3.0.

\code{details$savecall} determines whether the full call to \code{secr.fit} is 
saved in the output object. The default is TRUE except when called by 
\code{list.secr.fit} as names in the call are then evaluated, causing the 
output to become unwieldy.

\code{details$saveprogress} is a non-negative integer indicating the frequency
with which progress is saved to the RDS file named in details$progressfilename 
(overwritten each time).
Logical values FALSE and TRUE correspond to 0 (never) and 1 (after each likelihood 
evaluation). Less frequent saving is advised, as the operation is slow. If fitting
is terminated early it may be resumed from the saved file with function 
\code{\link{secr.refit}}. The value saved to the RDS file is a list of inputs 
for 'secr.fit'. The list has attribute 'log' that is a dataframe in which each 
row is a vector with the iteration number, log likelihood, time, and current 
values of the coefficients.

\code{details$progressfilename} is the name of an RDS file for 'saveprogress' 
(see previous). Ignored when saveprogress = 0. Default "progress.RDS".

\code{details$splitmarked} determines whether the home range centre of marked 
animals is allowed to move between the marking and sighting phases of a spatial 
capture--mark--resight study. The default is to assume a common home-range centre 
(\code{splitmarked = FALSE}). 

\code{details$telemetrytype} determines how telemetry data in the
attribute `xylist' are treated. `none' causes the xylist data to be
ignored. `dependent' uses information on the sampling distribution of
each home-range centre in the SECR likelihood. `concurrent' does that
and more: it splits capthist according to telemetry status and appends
all-zero histories to the telemetry part for any animals present in
xylist. The default is `concurrent'. 

\code{details$usecov} selects the mask covariate to be used for
normalization. NULL limits denominator for normalization to
distinguishing habitat from non-habitat.

\code{details$userDfn} is a user-provided function for modelling a density 
surface. See \href{https://www.otago.ac.nz/density/pdfs/secr-densitysurfaces.pdf}{secr-densitysurfaces.pdf}

\code{details$userdist} is either a function to compute non-Euclidean
distances between detectors and mask points, or a pre-computed matrix of
such distances. The first two arguments of the function should be
2-column matrices of x-y coordinates (respectively \eqn{k} detectors and
\eqn{m} mask points). The third argument is a habitat mask that defines
a non-Euclidean habitat geometry (a linear geometry is described in
documentation for the package `secrlinear'). The matrix
returned by the function must have exactly \eqn{k} rows and \eqn{m}
columns. When called with no arguments the function should return a
character vector of names for the required covariates of `mask',
possibly including the dynamically computed density `D` and a parameter
`noneuc' that will be fitted. A slightly expanded account is at
\link{userdist}, and full documentation is in the separate
document \href{https://www.otago.ac.nz/density/pdfs/secr-noneuclidean.pdf}{secr-noneuclidean.pdf}.

**Do not use `userdist' for polygon or transect detectors**

}

\references{

Efford, M. G., Dawson, D. K., Jhala, Y. V. and Qureshi, Q. (2016) 
Density-dependent home-range size revealed by spatially explicit 
capture--recapture. \emph{Ecography} \bold{39}, 676--688.

Efford, M. G. and Mowat, G. (2014) Compensatory heterogeneity in
capture--recapture data.\emph{Ecology} \bold{95}, 1341--1348.  

Gardner, B., Royle, J. A. and Wegan, M. T. (2009) Hierarchical models
for estimating density from DNA mark-recapture studies. \emph{Ecology}
\bold{90}, 1106--1115.

Royle, J. A., Chandler, R. B., Sun, C. C. and Fuller, A. K. (2013)
Integrating resource selection information with spatial
capture--recapture. \emph{Methods in Ecology and Evolution} \bold{4}, 520--530.

}

\seealso{
  \code{\link{secr.fit} }, \link{userdist}
}

\examples{

\dontrun{

## Demo of miscparm and userdist
## We fix the usual 'sigma' parameter and estimate the same 
## quantity as miscparm[1]. Differences in CI reflect the implied use 
## of the identity link for miscparm[1]. 

mydistfn3 <- function (xy1,xy2, mask) {
    if (missing(xy1)) return(character(0))
    xy1 <- as.matrix(xy1)
    xy2 <- as.matrix(xy2)
    miscparm <- attr(mask, 'miscparm')
    distmat <- edist(xy1,xy2) / miscparm[1]
    distmat
}

fit0 <- secr.fit (captdata)
fit <- secr.fit (captdata, fixed = list(sigma=1), details = 
    list(miscparm = c(sig = 20), userdist = mydistfn3))    
predict(fit0)
coef(fit)

}

}
