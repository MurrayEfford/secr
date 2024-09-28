\name{score.test}
\alias{score.test}
\alias{score.table}

\title{ Score Test for SECR Models }
\description{
  Compute score tests comparing a fitted model and a more general alternative model.
}
\usage{
score.test(secr, ..., betaindex = NULL, trace = FALSE, ncores = NULL, .relStep = 0.001,
    minAbsPar = 0.1)

score.table(object, ..., sort = TRUE, dmax = 10)

}

\arguments{
  \item{secr}{ fitted secr model }
  \item{...}{ one or more alternative models OR a fitted secr model }
  \item{trace}{ logical. If TRUE then output one-line summary at each
    evaluation of the likelihood }
  \item{ncores}{integer number of threads for parallel processing} 
  \item{.relStep}{see \code{\link[nlme]{fdHess}}}
  \item{ minAbsPar}{see \code{\link[nlme]{fdHess}}}
  \item{betaindex}{ vector of indices mapping fitted values to parameters in the alternative model }

  \item{object}{ score.test object or list of such objects }
  \item{sort}{ logical for whether output rows should be in descending order of AIC }
  \item{dmax}{ threshold of dAIC for inclusion in model set }

}
\details{

Score tests allow fast model selection (e.g. Catchpole & Morgan 1996).
Only the simpler model need be fitted. This implementation uses the
observed information matrix, which may sometimes mislead (Morgan et al.
2007). The gradient and second derivative of the likelihood function are
evaluated numerically at the point in the parameter space of the second
model corresponding to the fit of the first model. This operation uses
the function \code{fdHess} of the \pkg{nlme} package; the likelihood
must be evaluated several times, but many fewer times than would be
needed to fit the model. The score statistic is an approximation to the
likelihood ratio; this allows the difference in AIC to be estimated.

Covariates are inferred from components of the reference model
\code{secr}. If the new models require additional covariates these may
usually be added to the respective component of \code{secr}. 

Mapping of parameters between the fitted and alternative models
sometimes requires user intervention via the \code{betaindex} argument.
For example \code{betaindex} = c(1,2,4) is the correct mapping when
comparing the null model (D\eqn{\sim{~}}{~}1, g0\eqn{\sim{~}}{~}1,
sigma\eqn{\sim{~}}{~}1) to one with a behavioural effect on g0
(D\eqn{\sim{~}}{~}1, g0\eqn{\sim{~}}{~}b, sigma\eqn{\sim{~}}{~}1).

The arguments \code{.relStep} and \code{minAbsPar} control the numerical
gradient calculation and are passed directly to
\code{\link[nlme]{fdHess}}. More investigation is needed to determine
optimal settings.

\code{score.table} summarises one or more score tests in the form of a
model comparison table. The \dots argument here allows the inclusion of
additional score test objects (note the meaning differs from
\code{score.test}). Approximate AIC values are used to compute relative
AIC model weights for all models within dmax AIC units of the best
model.

If \code{ncores = NULL} then the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS is used (see \code{\link{setNumThreads}}).

}
\value{
  An object of class `score.test' that inherits from `htest', a list with components
  \item{statistic }{the value the chi-squared test statistic (score statistic)}
  \item{parameter }{degrees of freedom of the approximate chi-squared distribution of the test statistic (difference in number of parameters H0, H1)}
  \item{p.value }{probability of test statistic assuming chi-square distribution}
  \item{method }{a character string indicating the type of test performed }
  \item{data.name }{character string with null hypothesis, alternative hypothesis and arguments to function call from fit of H0 }
  \item{H0 }{simpler model }
  \item{np0 }{number of parameters in simpler model }
  \item{H1 }{alternative model }
  \item{H1.beta }{coefficients of alternative model }
  \item{AIC }{Akaike's information criterion, approximated from score statistic }
  \item{AICc }{AIC with small-sample adjustment of Hurvich & Tsai 1989 }

  If \dots defines several alternative models then a list of score.test objects is returned.

  The output from \code{score.table} is a dataframe with one row per model, including the reference model. 

}

\note{

This implementation is experimental. The AIC values, and values derived
from them, are approximations that may differ considerably from AIC
values obtained by fitting and comparing the respective models. Use of
the observed information matrix may not be optimal.

Weights were based on AICc rather than AIC prior to version 5.0.0.

}

\references{ 
Catchpole, E. A. and Morgan, B. J. T. (1996) Model selection of ring-recovery models using score tests. \emph{Biometrics} \bold{52}, 664--672.

Hurvich, C. M. and Tsai, C. L. (1989) Regression and time series model selection in small samples. \emph{Biometrika} \bold{76}, 297--307.

McCrea, R. S. and Morgan, B. J. T. (2011) Multistate mark-recapture
  model selection using score tests. \emph{Biometrics} \bold{67},
  234--241.
  
Morgan, B. J. T., Palmer, K. J. and Ridout, M. S. (2007) Negative score test statistic. \emph{American statistician} \bold{61}, 285--288.
}

\seealso{ 

\code{\link{AIC}}, 
\code{\link{LR.test}} 

}

\examples{

\dontrun{
    AIC (secrdemo.0, secrdemo.b)
    st <- score.test (secrdemo.0, g0 ~ b)
    st
    score.table(st)

    ## adding a time covariate to separate occasions (1,2) from (3,4,5)
    secrdemo.0$timecov <- data.frame(t2 = factor(c(1,1,2,2,2)))
    st2 <- score.test (secrdemo.0, g0 ~ t2)
    score.table(st,st2)
}

}
