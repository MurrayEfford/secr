\name{secr.fit}
\alias{secr.fit}

\title{ Spatially Explicit Capture--Recapture }

\description{

Estimate animal population density with data from an array of passive
detectors (traps) by fitting a spatial detection model by maximizing the
likelihood. Data must have been assembled as an object of class
\code{capthist}. Integration is by summation over the grid of points in
\code{mask}.

}

\usage{

secr.fit (capthist, model = list(D~1, g0~1, sigma~1), mask = NULL, buffer = NULL,
    CL = FALSE, detectfn = NULL, binomN = NULL, start = NULL, link = list(),
    fixed = list(), timecov = NULL, sessioncov = NULL, hcov = NULL,
    groups = NULL, dframe = NULL, details = list(), method =
    "Newton-Raphson", verify = TRUE, biasLimit = 0.01, trace = NULL,
    ncores = NULL, ...)

}

\arguments{
  \item{capthist}{ \code{\link{capthist}} object including capture data and detector (trap) layout }
  \item{mask}{ \code{\link{mask}} object or (for a multi-session analysis) a list of \code{mask} objects, one for each session}
  \item{buffer}{ scalar mask buffer radius if \code{mask} not specified
  (default 100 m)}
  \item{CL}{ logical, if true then the model is fitted by maximizing the conditional likelihood }
  \item{detectfn}{integer code or character string for shape of
  detection function 0 = halfnormal, 1 = hazard rate etc. -- see
  \link{detectfn}}
  \item{binomN}{ integer code for distribution of counts (see Details) }
  \item{start}{ vector of initial values for beta parameters, or \code{secr} object from which they may be derived }
  \item{link}{ list with optional components corresponding to `real'
  parameters (e.g., `D', `g0', `sigma'), each a character string in
  \{"log", "logit", "identity", "sin"\} for the link function of one real parameter }
  \item{fixed}{ list with optional components corresponding to real parameters giving the scalar value to which the parameter is to be fixed }
  \item{model}{ list with optional components each symbolically defining a linear predictor for one real parameter using \code{formula} notation }
  \item{timecov}{ optional dataframe of values of time (occasion-specific) covariate(s). }
  \item{sessioncov}{ optional dataframe of values of session-specific covariate(s). }
  \item{hcov}{ character name of individual covariate for known
  membership of mixture classes. }
  \item{groups}{ optional vector of one or more variables with which to form groups. Each element should be the name of a factor variable in the \code{covariates} attribute of \code{capthist}. }
  \item{dframe}{ optional data frame of design data for detection parameters }
  \item{details}{ list of additional settings, mostly model-specific (see Details) }
  \item{method}{ character string giving method for maximizing log likelihood }
  \item{verify}{ logical, if TRUE the input data are checked with \code{\link{verify}} }
  \item{biasLimit}{ numeric threshold for predicted relative bias due to
  buffer being too small }
  \item{trace}{ logical, if TRUE then output each evaluation of the likelihood, and other messages}
  \item{ncores}{ integer number of threads to use for parallel processing}
  \item{\dots}{ other arguments passed to the maximization function }
}

\details{

\code{secr.fit} fits a SECR model by maximizing the likelihood. The
likelihood depends on the detector type ("multi", "proximity", "count",
"polygon" etc.) of the \code{traps} attribute of \code{capthist}
(Borchers and Efford 2008, Efford, Borchers and Byrom 2009, Efford,
Dawson and Borchers 2009, Efford 2011). The `multi' form of the
likelihood is also used, with a warning, when detector type = "single"
(see Efford et al. 2009 for justification).

The default \code{model} is null (\code{model = list(D~1, g0~1,
sigma~1)} for \code{detectfn = 'HN'} and \code{CL = FALSE}), meaning
constant density and detection probability). The set of variables
available for use in linear predictors includes some that are
constructed automatically (t, T, b, B, bk, Bk, k, K), group (g), and
others that appear in the \code{covariates} of the input data. See also
\code{\link{usage}} for varying effort, \code{\link{timevaryingcov}} to
construct other time-varying detector covariates, and \href{https://www.otago.ac.nz/density/pdfs/secr-models.pdf}{secr-models.pdf}
and \href{../doc/secr-overview.pdf}{secr-overview.pdf} for more on
defining models.

\code{buffer} and \code{mask} are alternative ways to define the region
of integration (see \link{mask}). If \code{mask} is not specified then a
mask of type "trapbuffer" will be constructed automatically using the
specified buffer width in metres.

\code{hcov} is used to define a hybrid mixture model, used especially to
model sex differences (see \code{\link{hcov}}). (Allows some animals to
be of unknown class).

The length of \code{timecov} should equal the number of sampling
occasions (\code{ncol(capthist)}). Arguments \code{timecov},
\code{sessioncov} and \code{groups} are used only when needed for terms
in one of the model specifications. Default \code{link} is \code{list(D="log",
g0="logit", sigma="log")}.

If \code{start} is missing then \code{\link{autoini}} is used for D, g0
and sigma, and other beta parameters are set initially to arbitrary
values, mostly zero. \code{start} may be a previously fitted model. In
this case, a vector of starting beta values is constructed from the old
(usually nested) model and additional betas are set to zero. Mapping of
parameters follows the default in \code{\link{score.test}}, but user
intervention is not allowed. From 2.10.0 the new and old models need not
share all the same `real' parameters, but any new real parameters, such
as `pmix' for finite mixture models, receive a starting value of 0 on
the link scale (remembering e.g., invlogit(0) = 0.5 for parameter `pmix').

\code{binomN} (previously a component of \code{details}) determines the
distribution that is fitted for the number of detections of an individual
at a particular detector, on a particular occasion, when the detectors
are of type `count', `polygon' or `transect':

\itemize{
\item{binomN > 1  --- binomial with size binomN}
\item{binomN = 1  --- binomial with size determined by \link{usage}}
\item{binomN = 0  --- Poisson}
%\item{binomN < 0} {negative binomial with size abs(binomN) -- see 
%  \code{\link[stats]{NegBinomial}}}
}

The default with these detectors is to fit a Poisson distribution.
%The `size' parameter of the negative binomial is not estimated: it must be
%supplied. \code{binomN} should be an integer unless negative.

\code{details} is used for various specialized settings listed below. These are
described separately - see \code{\link{details}}.
\tabular{ll}{
autoini \tab session to use for starting values (default 1) \cr
centred \tab centre x-y coordinates\cr
chat \tab overdispersion of sighting counts Tu, Tm \cr
chatonly \tab compute overdispersion for Tu and Tm, then exit \cr
contrasts \tab coding of factor predictors \cr
convexpolygon \tab allows non-convex polygons (slower) \cr
Dfn \tab reparameterization of density model (seldom used directly) \cr
Dlambda \tab switch density reparameterization to trend model \cr
distribution \tab binomial vs Poisson N\cr
fastproximity \tab special handling of binary proximity detectors \cr
fixedbeta \tab specify fixed beta parameter(s)\cr
grain \tab grain argument of RcppParallel::parallelFor \cr
hessian \tab variance method\cr
ignoreusage \tab override usage in traps object of capthist\cr
intwidth2 \tab controls optimise when only one parameter\cr
knownmarks \tab known or unknown number of marked animals in sighting-only model \cr
LLonly \tab compute one likelihood for values in \code{start} \cr
maxdistance \tab distance threshold for selective mask \cr
miscparm \tab starting values for extra parameters fitted via userdist function \cr
newdetector \tab detector type to override detector(traps(capthist)) \cr
nsim \tab number of simulations to compute overdispersion \cr
param \tab optional parameterisation code\cr
relativeD \tab optional relative density conditional on \eqn{n} \cr
savecall \tab optionally suppress saving of call\cr
telemetrytype \tab treat telemetry data as independent, dependent or concurrent\cr
normalize \tab rescale detection to individual range use \cr
usecov \tab spatial covariate of use for normalization \cr
userdist \tab user-provided distance function or matrix \cr
%maxcallsize \tab limit to size of saved call \cr
}

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

A mark-resight model is fitted if the \code{\link{markocc}} attribute of the capthist 
`traps' object includes sighting occasions. See the vignette 
\href{https://www.otago.ac.nz/density/pdfs/secr-markresight.pdf}{secr-markresight.pdf} 
for a full account.

If \code{method = "Newton-Raphson"} then \code{\link[stats]{nlm}} is
used to maximize the log likelihood (minimize the negative log
likelihood); otherwise \code{\link[stats]{optim}} is used with the
chosen method ("BFGS", "Nelder-Mead", etc.).  If maximization fails a
warning is given appropriate to the method.

From secr 2.5.1, \code{method = "none"} may be used to skip likelihood
maximization and compute only the hessian for the current dataset at the
values in start, and the corresponding variance-covariance matrix of
beta parameters. The computation uses fdHess from \pkg{nlme}.

If \code{verify} = TRUE then \code{\link{verify}} is called to check
capthist and mask; analysis is aborted if "errors" are found. Some
conditions that trigger an "error" are benign (e.g., no detections in
some sessions of a multi-session study of a sparse population); use
\code{verify = FALSE} to avoid the check. See also Note.

If \code{buffer} is used rather than \code{mask}, and \code{biasLimit}
    is valid, then the estimated density is checked for bias due to the
    choice of buffer. A warning is generated when \code{buffer} appears
    to be too small (predicted RB(D-hat) > \code{biasLimit}, default 1\%
    relative bias). The prediction uses \code{\link{bias.D}}. No check
    is performed when \code{mask} is specified, when \code{biasLimit} is 0,
    negative or NA, or when the detector type is "polygon", "transect",
    "polygonX" or "transectX".

Function \code{\link{list.secr.fit}} is a way to fit several models at once.

}

\value{

When \code{details$LLonly = TRUE} a single log-likelihood is returned, with attributes 

\item{npar}{ number of parameters to be estimated,} 
\item{preptime}{ elapsed setup time in seconds,}
\item{LLtime}{ elapsed time for single likelihood evaluation, exclusive of setup.}

Otherwise, \code{secr.fit} returns an object of class secr representing the fitted SECR model.  This has
components

  \item{call }{function call) }

  \item{capthist }{saved input}
  \item{mask }{saved input}
  \item{detectfn }{saved input}
  \item{CL }{saved input}
  \item{timecov }{saved input}
  \item{sessioncov }{saved input}
  \item{hcov }{saved input}
  \item{groups }{saved input}
  \item{dframe }{saved input}
  \item{designD }{design matrix for density model; may be NULL}
  \item{designNE }{design matrix for noneuc model; may be NULL}
  \item{design }{reduced design matrices for detection parameters, parameter table and parameter
    index array for actual animals (see \code{\link{secr.design.MS}})}
  \item{design0 }{reduced design matrices for detection parameters, parameter table and parameter
    index array for `naive' animal (see \code{\link{secr.design.MS}})}

  \item{start }{vector of starting values for beta parameters} 

  \item{link }{list with one component for each real parameter
  (typically `D', `g0', `sigma'),giving the name of the link function
  used for each real parameter.}

  \item{fixed }{saved input}  

  \item{parindx }{list with one component for each real parameter giving
    the indices of the `beta' parameters associated with each real
    parameter}  

  \item{model }{saved input}

  \item{details }{saved input}

  \item{vars }{vector of unique variable names in \code{model} }

  \item{betanames }{names of beta parameters}

  \item{realnames }{names of fitted (real) parameters }

  \item{fit }{list describing the fit (output from \code{nlm} or
  \code{optim}) }

  \item{beta.vcv }{variance-covariance matrix of beta parameters }  
  
  \item{smoothsetup}{list of objects specifying smooths in \pkg{mgcv}}
  \item{learnedresponse}{ logical; TRUE if any learned response in detection model}

  \item{version }{secr version number }
  \item{starttime }{character string of date and time at start of fit }
  \item{proctime }{processor time for model fit, in seconds }
  
The environment variable RCPP_PARALLEL_NUM_THREADS is updated if an integer value is provided for \code{ncores}.

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area searches. \emph{Ecology}
\bold{92}, 2202--2207.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture--recapture: likelihood-based methods. In:
D. L. Thompson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer. Pp. 255--269.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\note{

One system of units is used throughout \pkg{secr}. Distances are in metres and
areas are in hectares (ha). The unit of density is animals per
hectare. 1 ha = 10000 m^2 = 0.01 km^2. To convert density to animals /
km^2, multiply by 100.

When you display an `secr' object by typing its name at the command
prompt, you implicitly call its `print' method \code{\link{print.secr}}, which
in turn calls \code{\link{predict.secr}} to tabulate estimates of the `real'
parameters. Confidence limits (lcl, ucl) are for a 100(1-alpha)\%
interval, where alpha defaults to 0.05 (95\% interval); alpha may be
varied in \code{print.secr} or \code{predict.secr}.

\code{AIC}, \code{logLik} and \code{vcov} methods are also
provided. Take care with using AIC: not all models are comparable (see
Notes section of \code{\link{AIC.secr}}) and large differences in AIC
may relate to trivial differences in estimated density.

\code{derived} is used to compute the derived parameters `esa'
(effective sampling area) and `D' (density) for models fitted by
maximizing the conditional likelihood (CL = TRUE).

Components `version' and `starttime' were introduced in version 1.2.7,
and recording of the completion time in `fitted' was discontinued.

 The Newton-Raphson algorithm is fast, but it sometimes fails to compute
the information matrix correctly, causing some or all standard errors to
be set to NA. This usually indicates a major problem in fitting the
model, and parameter estimates should not be trusted. See
\link{Troubleshooting}.

The component D in output was replaced with N from version 2.3. Use
\code{\link{region.N}} to obtain SE or confidence intervals for N-hat,
or to infer N for a different region.

Prior to version 2.3.2 the buffer bias check could be switched off by
setting \code{verify = FALSE}. This is now done by setting
\code{biasLimit = 0} or \code{biasLimit = NA} .

}

\seealso{

  \link{Detection functions},
  \code{\link{AIC.secr}},
  \code{\link{capthist}},
  \code{\link{details}},
  \code{\link{derived}},
  \code{\link{hcov}},
  \code{\link{mask}},
  \code{\link{list.secr.fit}},
  \code{\link{predict.secr}},
  \code{\link{print.secr}},
  \code{\link{region.N}},
  \link{Speed tips}
  \link{Troubleshooting}
  \link{userdist}
  \code{\link{usage}},
  \code{\link{vcov.secr}},
  \code{\link{verify}},

}

\section{Warning}{** Mark-resight data formats and models are experimental in secr 2.10.0 and subject to change **}

\examples{

\dontrun{

## construct test data (array of 48 `multi-catch' traps)

detectors <- make.grid (nx = 6, ny = 8, detector = "multi")
detections <- sim.capthist (detectors, popn = list(D = 10,
    buffer = 100), detectpar = list(g0 = 0.2, sigma = 25))

## fit & print null (constant parameter) model
secr0 <- secr.fit (detections)
secr0   ## uses print method for secr

## compare fit of null model with learned-response model for g0

secrb <- secr.fit (detections, model = g0~b)
AIC (secr0, secrb)

## typical result

##                  model   detectfn npar    logLik     AIC    AICc dAICc  AICwt
## secr0 D~1 g0~1 sigma~1 halfnormal    3 -347.1210 700.242 700.928 0.000 0.7733
## secrb D~1 g0~b sigma~1 halfnormal    4 -347.1026 702.205 703.382 2.454 0.2267
}
}

\keyword{ models }

