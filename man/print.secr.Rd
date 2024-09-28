\name{print.secr}
\alias{print.secr}
\alias{summary.secr}

\title{ Print or Summarise secr Object }

\description{

  Print results from fitting a spatially explicit capture--recapture model or
  generate a list of summary values.

}

\usage{
\method{print}{secr} (x, newdata = NULL, alpha = 0.05, deriv = FALSE, call = TRUE, ...)
\method{summary}{secr} (object, newdata = NULL, alpha = 0.05, deriv = FALSE, ...)
}

\arguments{
  \item{x}{ \code{secr} object output from \code{secr.fit} }
  \item{object}{ \code{secr} object output from \code{secr.fit} }
  \item{newdata}{ optional dataframe of values at which to evaluate model }
  \item{alpha}{ alpha level }
  \item{deriv}{ logical for calculation of derived D and esa }
  \item{call}{ logical; if TRUE the call is printed }
  \item{...}{ other arguments optionally passed to derived.secr }
}

\details{

Results from \code{print.secr} are potentially complex and depend upon the analysis (see
below). Optional \code{newdata} should be a dataframe with a column for
each of the variables in the model. If \code{newdata} is missing then a
dataframe is constructed automatically.  Default \code{newdata} are for
a naive animal on the first occasion; numeric covariates are set to zero
and factor covariates to their base (first) level. Confidence intervals
are 100 (1 -- alpha) \% intervals.

\tabular{ll}{
call            \tab the function call (optional) \cr
version,time    \tab secr version, date and time fitting started, and elapsed time \cr

Detector type   \tab `single', `multi', `proximity' etc. \cr
Detector number \tab number of detectors \cr
Average spacing \tab \cr
x-range         \tab \cr
y-range         \tab \cr
New detector type \tab as fitted when details$newdetector specified \cr
N animals       \tab number of distinct animals detected \cr
N detections    \tab number of detections \cr
N occasions     \tab number of sampling occasions \cr
Mask area       \tab \cr

Model          \tab model formula for each `real' parameter \cr
Fixed (real)   \tab fixed real parameters \cr
Detection fn   \tab detection function type (halfnormal or hazard-rate) \cr
N parameters   \tab number of parameters estimated  \cr
Log likelihood \tab log likelihood \cr
AIC            \tab Akaike's information criterion \cr
AICc           \tab AIC with small sample adjustment (Burnham and Anderson 2002) \cr
Beta parameters \tab coef of the fitted model, SE and confidence intervals \cr
vcov            \tab variance-covariance matrix of beta parameters  \cr
Real parameters \tab fitted (real) parameters evaluated at base levels of covariates \cr
Derived parameters \tab derived estimates of density and mean effective sampling area (optional) \cr
}

Derived parameters (see \code{\link{derived}}) are computed only if
\code{deriv = TRUE}.

}

\value{

The \code{summary} method constructs a list of outputs similar to those printed by the \code{print} method, but somewhat more concise and re-usable:

\tabular{ll}{
versiontime     \tab secr version, and date and time fitting started \cr

traps           \tab detector summary \cr
capthist        \tab capthist summary \cr
mask            \tab mask summary \cr
modeldetails    \tab miscellaneous model characteristics (CL etc.) \cr
AICtable        \tab single-line output of AIC.secr \cr
coef            \tab table of fitted coefficients with CI \cr
predicted       \tab predicted values (`real' parameter estimates) \cr
derived         \tab output of derived.secr (optional) \cr
}

}

\references{ 

Burnham, K. P. and Anderson, D. R. (2002) \emph{Model selection and multimodel inference: a practical information-theoretic approach}. Second edition. New York: Springer-Verlag.

}

\seealso{ \code{\link{AIC.secr}}, \code{\link{secr.fit}} }

\examples{

## load & print previously fitted null (constant parameter) model
print(secrdemo.0)

summary(secrdemo.0)

## combine AIC tables from list of summaries
do.call(AIC, lapply(list(secrdemo.b, secrdemo.0), summary))

\dontrun{

print(secrdemo.CL, deriv = TRUE)

}

}

