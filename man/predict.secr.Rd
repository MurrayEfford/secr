\name{predict.secr}
\alias{predict.secr}
\alias{predict.secrlist}
\alias{detectpar}
\alias{detectpar.secr}

\title{ SECR Model Predictions }

\description{
Evaluate a spatially explicit capture--recapture model. That is, compute the `real' parameters corresponding to the `beta' parameters of a fitted model for arbitrary levels of any variables in the linear predictor.
}

\usage{

\method{predict}{secr}(object, newdata = NULL, realnames = NULL, type = c("response", "link"), 
    se.fit = TRUE, alpha = 0.05, savenew = FALSE, \dots)

\method{predict}{secrlist}(object, newdata = NULL, realnames = NULL, type = c("response", "link"),
    se.fit = TRUE, alpha = 0.05, savenew = FALSE, \dots)

\method{detectpar}{secr}(object, \dots, byclass = FALSE) 

}

\arguments{
  \item{object}{ \code{secr} object output from \code{secr.fit}, or list
  of \code{secr} objects (secrlist) }
  \item{newdata}{ optional dataframe of values at which to evaluate model }
  \item{realnames}{ character vector of real parameter names }
  \item{type}{ character; type of prediction required. The default ("response") provides estimates of the `real' parameters.} 
  \item{se.fit}{ logical for whether output should include SE and confidence intervals }
  \item{alpha}{ alpha level for confidence intervals }
  \item{savenew}{ logical for whether newdata should be saved }
  \item{\dots}{ other arguments passed to \code{\link{makeNewData}} }
  \item{byclass}{ logical; if TRUE values are returned for each latent
  class in a mixture model, or class in a hybrid mixture (hcov) model}
}
\details{
  
The variables in the various linear predictors are described in
\href{https://www.otago.ac.nz/density/pdfs/secr-models.pdf}{secr-models.pdf} and listed for the particular model in the
\code{vars} component of \code{object}.

Optional \code{newdata} should be a dataframe with a column for each of
the variables in the model (see `vars' component of \code{object}).  If
\code{newdata} is missing then a dataframe is constructed automatically.

Default \code{newdata} are for a naive animal on the first occasion;
numeric covariates are set to zero and factor covariates to their base
(first) level. From secr 3.1.4 the argument `all.levels' may be passed 
to \code{\link{makeNewData}}; if TRUE then the default newdata includes 
all factor levels.

\code{realnames} may be used to select a subset of parameters.

Standard errors for parameters on the response (real) scale are by the
delta method (Lebreton et al. 1992), and confidence intervals are
backtransformed from the link scale. 

The value of \code{newdata} is optionally saved as an attribute.

\code{detectpar} is used to extract the detection parameter estimates
from a simple model to pass to functions such as
\code{esa.plot}. \code{detectpar} calls \code{predict.secr}. Parameters
will be evaluated by default at base levels of the covariates, although
this may be overcome by passing a one-line \code{newdata} to
\code{predict} via the \dots argument. Groups and mixtures are a
headache for \code{detectpar}: it merely returns the estimated detection
parameters of the first group or mixture.

If the `a0' parameterization has been used in \code{secr.fit} (i.e.,
\code{object$details$param == 3}) then \code{detectpar} automatically
backtransforms (a0, sigma) to (g0, sigma) or (lambda0, sigma) depending
on the value of \code{object$detectfn}.

}

\value{
  When \code{se.fit} = FALSE, a dataframe identical to \code{newdata} except for the addition of one column for each `real' parameter. Otherwise, a list with one component for each row in \code{newdata}. Each component is a dataframe with one row for each `real' parameter (density, g0, sigma, b) and columns as below

\tabular{ll}{
link        \tab link function \cr 
estimate    \tab estimate of real parameter \cr
SE.estimate \tab standard error of the estimate \cr
lcl         \tab lower 100(1--alpha)\% confidence limit \cr
ucl         \tab upper 100(1--alpha)\% confidence limit \cr
}

When \code{newdata} has only one row, the structure of the list is
`dissolved' and the return value is one data frame.

For \code{detectpar}, a list with the estimated values of detection
parameters (e.g., g0 and sigma if detectfn = "halfnormal"). In the case
of multi-session data the result is a list of lists (one list per
session).

}

\note{

  \code{\link{predictDsurface}} should be used for predicting density at many
  points from a model with spatial variation. This deals automatically
  with scaling of x- and y-coordinates, and is much is faster than
  predict.secr. The resulting Dsurface object has its own plot method.
  
  The argument `scaled' was removed from both predict methods in version 2.10 as the scaleg0 and scalesigma features had been superceded by other parameterisations.
  
  Overdispersion results in confidence intervals that are too narrow. See \code{\link{adjustVarD}} for a partial solution.
  
}

\references{ 
Lebreton, J.-D., Burnham, K. P., Clobert, J. and Anderson, D. R. (1992) Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. \emph{Ecological Monographs} \bold{62}, 67--118. }

\seealso{ 

\code{\link{secr.fit}}, 
\code{\link{predictDsurface}},  
\code{\link{adjustVarD}}, 
\code{\link{makeNewData}}

}

\examples{

## load previously fitted secr model with trap response
## and extract estimates of `real' parameters for both
## naive (b = 0) and previously captured (b = 1) animals

predict (secrdemo.b, newdata = data.frame(b = 0:1))

## OR from secr 3.1.4 
predict (secrdemo.b, all.levels = TRUE)

temp <- predict (secrdemo.b, all.levels = TRUE, save = TRUE)
attr(temp, "newdata")

detectpar(secrdemo.0)
}

\keyword{ models }

