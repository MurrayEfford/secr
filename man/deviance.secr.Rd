\name{deviance}
\alias{deviance}
\alias{df.residual}
\alias{deviance.secr}
\alias{df.residual.secr}

\title{
Deviance of fitted secr model and residual degrees of freedom
}

\description{

Compute the deviance or residual degrees of freedom of a fitted secr
model, treating multiple sessions and groups as independent. The
likelihood of the saturated model depends on whether the `conditional'
or `full' form was used, and on the distribution chosen for the number
of individuals observed (Poisson or binomial).

}
\usage{
\method{deviance}{secr}(object, ...)
\method{df.residual}{secr}(object, ...)
}

\arguments{
  \item{object}{ secr object from secr.fit }
  \item{\dots}{ other arguments (not used)}
}
\details{

The deviance is \eqn{-2log(\hat{L}) + 2log(L_{sat})}{-2.log(L-hat) +
2.log(L_sat)}, where \eqn{\hat{L}}{L-hat} is the value of the
log-likelihood evaluated at its maximum, and \eqn{L_{sat}}{L_sat} is the
log-likelihood of the saturated model, calculated thus:

Likelihood conditional on \eqn{n} -

\eqn{L_{sat} = \log(n!) + \sum\limits _{\omega} [n_\omega \log (\frac{n_\omega}{n}) - \log 
 (n_\omega !)]}{ L_sat = log(n!) + sum-over-w ( n_w log (n_w / n) - log (n_w !)))}

Full likelihood, Poisson \eqn{n} -

\eqn{L_{sat} = n\log(n) - n + \sum\limits _{\omega} [n_\omega \log (\frac{n_\omega}{n}) - \log (n_\omega !)]}{ L_sat = nlog(n) - n + sum-over-w [n_w log (n_w / n) - log (n_w !)] }

Full likelihood, binomial \eqn{n} -

\eqn{L_{sat} = n\log(\frac{n}{N}) + (N-n)\log(\frac{N-n}{N}) + \log (\frac{N!}{(N-n)!}) + \sum\limits _{\omega} [n_\omega \log (\frac{n_\omega}{n}) - \log (n_\omega !)]}{ L_sat = nlog(n/N) + (N-n)log((N-n)/N) + log(N!/(N-n)!) + sum-over-w [n_w log (n_w/n) - log (n_w!)]}

\eqn{n} is the number of individuals observed at least once, \eqn{n_\omega}{n_w} is the number of distinct histories, and \eqn{N} is the number in a chosen area \eqn{A} that we estimate by \eqn{\hat{N} = \hat{D}A}{N-hat = D-hat.A}.

The residual degrees of freedom is the number of distinct detection
histories minus the number of parameters estimated. The detection
histories of two animals are always considered distinct if they belong to
different groups.

When samples are (very) large the deviance is expected to be distributed
as \eqn{\chi^2}{chi-square} with \eqn{n_\omega - p}{n_w - p} degrees of
freedom when \eqn{p} parameters are estimated. In reality, simulation is
needed to assess whether a given value of the deviance indicates a
satisfactory fit, or to estimate the overdispersion parameter \eqn{c}.
\code{sim.secr} is a convenient tool. 

}

\value{

The scalar numeric value of the deviance or the residual degress of freedom extracted
from the fitted model.

}
\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

}

\seealso{ 

\code{\link{secr.fit}}, \code{\link{sim.secr}}

}

\examples{

deviance(secrdemo.0)
df.residual(secrdemo.0)

}


