\name{confint.secr} 
\alias{confint.secr}

\title{ Profile Likelihood Confidence Intervals }

\description{

Compute profile likelihood confidence intervals for `beta' or `real' parameters of a
spatially explicit capture-recapture model,

}

\usage{

\method{confint}{secr} (object, parm, level = 0.95, newdata = NULL,
tracelevel = 1, tol = 0.0001, bounds = NULL, ncores = NULL, ...)

}

\arguments{
  \item{object}{ \code{secr} model object }
  \item{parm}{ numeric or character vector of parameters }
  \item{level}{ confidence level (1 -- alpha) }
  \item{newdata}{ optional dataframe of values at which to evaluate model }
  \item{tracelevel}{ integer for level of detail in reporting (0,1,2) }
  \item{tol}{ absolute tolerance (passed to \link{uniroot}) }
  \item{bounds}{ numeric vector of outer starting values -- optional }
  \item{ncores}{ number of threads used for parallel processing }
  \item{\dots}{ other arguments (not used) }
}
\details{

If \code{parm} is numeric its elements are interpreted as the indices of
`beta' parameters; character values are interpreted as `real'
parameters. Different methods are used for beta parameters and real
parameters. Limits for the \eqn{j}-th beta parameter are found by a
numerical search for the value satisfying \eqn{-2(l_j(\beta_j) - l) =
q}{-2(l_j(beta_j) - l) = q}, where \eqn{l} is the maximized log
likelihood, \eqn{l_j(\beta_j)}{l_j(beta_j)} is the maximized profile log
likelihood with \eqn{\beta_j}{beta_j} fixed, and \eqn{q} is the
\eqn{100(1-\alpha)}{100(1-alpha)} quantile of the
\eqn{\chi^2}{chi-square} distribution with one degree of freedom. Limits
for real parameters use the method of Lagrange multipliers (Fletcher and
Faddy 2007), except that limits for constant real parameters are
backtransformed from the limits for the relevant beta parameter.

If \code{bounds} is provided it should be a 2-vector or matrix of 2
columns and length(parm) rows.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

} 

\value{ 

A matrix with one row for each parameter in \code{parm}, and columns
giving the lower (lcl) and upper (ucl) 100*\code{level}% limits.

} 

\references{ 

Evans, M. A., Kim, H.-M. and O'Brien, T. E. (1996) An application of
profile-likelihood based confidence interval to capture--recapture
estimators. \emph{Journal of Agricultural, Biological and Experimental
Statistics} \bold{1}, 131--140.

Fletcher, D. and Faddy, M. (2007) Confidence intervals for expected
abundance of rare species. \emph{Journal of Agricultural, Biological and
Experimental Statistics} \bold{12}, 315--324.

Venzon, D. J. and Moolgavkar, S. H. (1988) A method for computing
profile-likelihood-based confidence intervals. \emph{Applied Statistics}
\bold{37}, 87--94.

} 

\note{

Calculation may take a long time, so probably you will do it
only after selecting a final model. 

The \R function \code{\link{uniroot}} is used to search for the roots of
\eqn{-2(l_j(\beta_j) - l) = q}{-2(l_j(beta_j) - l) = q} within a
suitable interval. The interval is anchored at one end by the MLE, and
at the other end by the MLE inflated by a small multiple of the
asymptotic standard error (1, 2, 4 or 8 SE are tried in turn, using the
smallest for which the interval includes a valid solution). 

A more efficient algorithm was proposed by Venzon and Moolgavkar (1988);
it has yet to be implemented in \pkg{secr}, but see \code{plkhci} in
the package \pkg{Bhat} for another \R implementation.

}

\examples{

\dontrun{

## Limits for the constant real parameter "D"
confint(secrdemo.0, "D")   

}

}


