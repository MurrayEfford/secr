\name{mask.check}
\alias{mask.check}

\title{
Mask Diagnostics
}

\description{
\code{mask.check} evaluates the effect of varying buffer width and
    mask spacing on either the likelihood or density estimates from
    secr.fit().
}

\usage{
mask.check(object, buffers = NULL, spacings = NULL, poly = NULL,
    LLonly = TRUE, realpar = NULL, session = 1, file = NULL,
    drop = "", tracelevel = 0, ...)
}

\arguments{
  \item{object}{object of class `capthist' or `secr'}
  \item{buffers}{vector of buffer widths}
  \item{spacings}{vector of mask spacings}
  \item{poly}{ matrix of two columns, the x- and  y-coordinates of a bounding polygon (optional) }
  \item{LLonly}{logical; if TRUE then only the log likelihood is
    computed}
  \item{realpar}{list of parameter values}
  \item{session}{vector of session indices (used if \code{object} spans multiple sessions)}
  \item{file}{name of output file (optional)}
  \item{drop}{character vector: names of fitted secr object to omit}
  \item{tracelevel}{integer for level of detail in reporting (0,1,2)}
  \item{\dots}{other arguments passed to secr.fit}
}
\details{

  Masks of varying buffer width and spacing are constructed with the
  `trapbuffer' method in \code{make.mask}, using the detector locations
  (`traps') from either a capthist object or a previous execution of
  \code{secr.fit}. Default values are provided for \code{buffers} and
  \code{spacings} if \code{object} is of class `secr' (respectively c(1,
  1.5, 2) and c(1, 0.75, 0.5) times the values in the existing
  mask). The default for \code{buffers} will not work if a detector is
  on the mask boundary, as the inferred buffer is then 0.

  Variation in the mask may be assessed for its effect on --
  \itemize{   
      \item the log-likelihood evaluated for given values of the
      parameters (\code{LLonly = TRUE})
      \item estimates from maximizing the likelihood with each mask
	(\code{LLonly = FALSE})      
      }      

  \code{realpar} should be a list with one named component for each real
  parameter (see Examples). It is relevant only if \code{LLonly =
  TRUE}. \code{realpar} may be omitted if \code{object} is of class
  `secr'; parameter values are then extracted from \code{object}.

  \code{session} should be an integer or character vector suitable for
  indexing sessions in \code{object}, or in \code{object$capthist} if
  \code{object} is a fitted model. Each session is considered
  separately; a model formula that refers to session or uses session
  covariates will cause an error.
  
  If \code{file} is specified then detailed
  results (including each model fit when \code{LLonly = FALSE}) are
  saved to an external .RData file. Loading this file creates or
  overwrites object(s) in the workspace: \code{mask.check.output} if
  \code{LLonly = TRUE}, otherwise \code{mask.check.output} and
  \code{mask.check.fit}. For multiple sessions these are replaced by
  lists with one component per session (\code{mask.check.outputs} and
  \code{mask.check.fits}). The \code{drop} argument is passed to
  \code{\link{trim}} and applied to each fitted model; use it to save
  space, at the risk of limiting further computation on the fitted
  models.

  \code{tracelevel>0} causes more verbose reporting of progress during
  execution.
  
  The \dots argument may be used to override existing settings in
  \code{object} - for example, a conditional likelihood fit (\code{CL =
  T}) may be selected even if the original model was fitted by
  maximizing the full likelihood.
  
}

\value{
  
  Array of log-likelihoods (\code{LLonly = TRUE}) or estimates
  (\code{LLonly = FALSE}) for each combination of \code{buffers} and
  \code{spacings}. The array has 3 dimensions if \code{LLonly = FALSE}
  and both \code{buffers} and \code{spacings} have multiple levels;
  otherwise it collapses to a matrix. Rows generally represent
  \code{buffers}, but rows represent \code{spacings} if a single buffer
  is specified.
  
}

\note{

  When \code{LLonly = TRUE} the functionality of \code{mask.check}
  resembles that of the `Tools | ML SECR log likelihood' menu option in
  Density 5. The help page in Density 5 for ML SECR 2-D integration (see
  index) may be helpful.

  Warning messages from \code{secr.fit} are suppressed. `capthist'
  data provided via the \code{object} argument are checked with
  \code{\link{verify.capthist}} if \code{tracelevel > 0}.
  
  The likelihood-only method is fast, but not definitive. It is
  reasonable to aim for stability in the third decimal place of the log
  likelihood. Slight additional variation in the likelihood may cause
  little change in the estimates; the only way to be sure is to check
  these by setting \code{LLonly = FALSE}.

  The performance of a mask depends on the detection function; be sure
  to set the \code{detectfn} argument appropriately. The hazard rate
  function has a fat tail that can be problematic.
  
  When provided with an `secr' \code{object}, \code{mask.check}
  constructs a default vector of buffer widths as multiples of the
  buffer used in \code{object} \emph{even though that value is not saved
  explicitly}. For this trick, detector locations in
  \code{traps(object$capthist)} are compared to the bounding box of
  \code{object$mask}; the base level of buffer width is the maximum
  possible within the bounding box.
  
}

\section{Warning}{
  
  \code{mask.check()} may fail if \code{object} is a fitted `secr' model
  and a data object named in the original call of \code{secr.fit()}
  (i.e. \code{object$call}) is no longer in the working environment
  (\code{secr.fit} arguments capthist, mask, verify & trace are
  exempt). Fix by any of (1) applying \code{mask.check} directly to the
  `capthist' object, specifying other arguments (\code{buffers},
  \code{spacings}, \code{realpar}) as needed, (2) re-fitting the model
  and running \code{mask.check} in the same environment, (3) specifying
  the offending argument(s) in \dots, or (4) re-creating the required
  data objects(s) in the working environment, possibly from saved inputs
  in \code{object} (e.g., \code{mytimecov <- myfit$timecov}).

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially explicit
capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand. \url{https://www.otago.ac.nz/density/}.

}

\seealso{
\code{\link{esaPlot}}, \code{\link{make.mask}}, \code{\link{secr.fit}} 
}

\examples{

% not run because (1) slow (2) writes output file
\dontrun{

## from a capthist object, specifying almost everything
mask.check (possumCH, spacings = c(20, 30), buffers =c(200, 300),
    realpar = list(g0 = 0.2, sigma = 50), CL = TRUE)

## from a fitted model, using defaults
mask.check (stoat.model.HN)
## LL did not change with varying buffer (rows) or spacing (cols):
##         78.125  58.59375   39.0625
## 1000 -144.0015 -144.0015 -144.0015
## 1500 -144.0017 -144.0017 -144.0017
## 2000 -144.0017 -144.0017 -144.0017

## fit new models for each combination of buffer & spacing,
## and save fitted models to a file
mask.check (stoat.model.HN, buffers = 1500, spacings =
    c(40,60,80), LLonly = FALSE, file = "test", CL = TRUE)

## look in more detail at the preceding fits
## restores objects `mask.check.output' and `mask.check.fit'
load("test.RData")  
lapply(mask.check.fit, predict)
lapply(mask.check.fit, derived)

## multi-session data
mask.check(ovenbird.model.1, session = c("2005","2009"))

## clipping mask
txtfilename <- system.file("extdata/possumarea.txt", package = "secr")
possumarea <- read.table(txtfilename, header = TRUE)
mask.check (possum.model.0, spacings = c(20, 30), buffers =
    c(200, 300), poly = possumarea, LLonly = FALSE,
    file = "temp", CL = TRUE)

## review fitted models
load ("temp.RData")
par(mfrow = c(2,2), mar = c(1,4,4,4))
for (i in 1:4) {
    plot(traps(mask.check.fit[[i]]$capthist), border = 300,
        gridlines = FALSE)
    plot(mask.check.fit[[i]]$mask, add = TRUE)
    lines(possumarea)
    text ( 2698618, 6078427, names(mask.check.fit)[i])
    box()
}
par(mfrow = c(1,1), mar = c(5,4,4,2) + 0.1)    ## defaults

}

}

\keyword{ manip }
