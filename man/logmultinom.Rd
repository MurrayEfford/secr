\name{logmultinom}
\alias{logmultinom}
\title{
Multinomial Coefficient of SECR Likelihood
}
\description{
Compute the constant multinomial component of the SECR log likelihood
}
\usage{
logmultinom(capthist, grp = NULL)
}
\arguments{
  \item{capthist}{ \code{\link{capthist}} object }
  \item{grp}{ factor defining group membership, or a list (see Details) }
}
\details{

  For a particular dataset and grouping, the multinomial coefficient is
  a constant; it does not depend on the parameters and may be ignored
  when maximizing the likelihood to obtain parameter
  estimates. Nevertheless, the log likelihood reported by
  \code{secr.fit} includes this component \emph{unless} the
  detector type is `signal', `polygon', `polygonX', `transect' or
  `transectX' (from 2.0.0).

  If \code{grp} is NULL then all animals are assumed to belong to one
  group. Otherwise, the length of \code{grp} should equal the number of rows of
  \code{capthist}.

  \code{grp} may also be any vector that can be coerced
  to a factor. If \code{capthist} is a multi-session capthist object
  then \code{grp} should be a list with one factor per session.

  If capture histories are not assigned to groups the value is the
  logarithm of
  \deqn{{{n}\choose{n_1, ..., n_C}} = {{n!} \over {n_1! n_2! ... n_C!}} }{see
  pdf manual} where \eqn{n} is the total number of
  capture histories and \eqn{n_1} ... \eqn{n_C} are the frequencies with
  which each of the \eqn{C} unique capture histories were observed.

  If capture histories are assigned to \eqn{G} groups the value is the
  logarithm of \deqn{{ \prod_{g=1}^{G} {{n_g!} \over {n_{g1}! n_{g2}! ... n_{gC_g}}!}}  }{see
  pdf manual} where \eqn{n_g} is the number of capture histories of
  group \eqn{g} and \eqn{n_{g1}} ... \eqn{n_{gC_g}} are the frequencies with
  which each of the \eqn{C_g} unique capture histories were observed for
  group \eqn{g}.

  For multi-session data, the value is the sum of the single-session
  values. Both session structure and group structure therefore affect
  the value computed. Users will seldom need this function.

}
\value{
  The numeric value of the log likelihood component.
}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture--recapture: likelihood-based methods. In:
D. L. Thompson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer. Pp. 255--269.

}

\seealso{\code{\link{stoatDNA}}}

\examples{

## no groups
logmultinom(stoatCH)

}
