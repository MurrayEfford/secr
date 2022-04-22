\name{autoini}
\alias{autoini}

\title{ Initial Parameter Values for SECR }

\description{
  Find plausible initial parameter values for \code{\link{secr.fit}}.  A
  simple SECR model is fitted by a fast ad hoc method.
}

\usage{
autoini(capthist, mask, detectfn = 0, thin = 0.2, tol = 0.001, 
    binomN = 1, adjustg0 = TRUE, adjustsigma = 1.2, ignoreusage = FALSE, 
    ncores = NULL)
}

\arguments{

\item{capthist}{\code{capthist} object}

\item{mask}{\code{mask} object compatible with the detector layout in
    \code{capthist}}

\item{detectfn}{integer code or character string for shape of
  detection function 0 = halfnormal}

\item{thin}{proportion of points to retain in mask }

\item{tol}{numeric absolute tolerance for numerical root finding }

\item{binomN}{integer code for distribution of counts (see \code{\link{secr.fit}}) }

\item{adjustg0}{logical for whether to adjust g0 for usage (effort) and
  binomN }

\item{adjustsigma}{numeric scalar applied to RPSV(capthist, CC = TRUE) }

\item{ignoreusage}{logical for whether to discard usage information from
\code{traps(capthist)}}

\item{ncores}{ integer number of threads to be used for parallel processing}

}

\details{

Plausible starting values are needed to avoid numerical
problems when fitting SECR models.  Actual models
to be fitted will usually have more than the three basic parameters
output by \code{autoini}; other initial values can usually be set to
zero for \code{secr.fit}.  If the algorithm encounters problems obtaining
a value for g0, the default value of 0.1 is returned.

Only the halfnormal detection function is currently available in \code{autoini} (cf
other options in e.g. \link{detectfn} and \code{\link{sim.capthist}}).

\code{autoini} implements a modified version of the algorithm proposed
by Efford et al. (2004).  In outline, the algorithm is
\enumerate{
\item Find value of sigma that predicts the 2-D dispersion of individual locations (see \code{\link{RPSV}}). 
\item Find value of g0 that, with sigma, predicts the observed mean number of captures per individual (by algorithm of Efford et al. (2009, Appendix 2))
\item Compute the effective sampling area from g0, sigma, using thinned mask (see \code{\link{esa}})
\item Compute D = \eqn{n}/esa(g0, sigma), where \eqn{n} is the number of individuals detected
}

Here `find' means solve numerically for zero difference between the observed and predicted values, using \code{\link{uniroot}}.

Halfnormal sigma is estimated with \code{RPSV(capthist, CC = TRUE)}.  The factor \code{adjustsigma} is applied as a crude correction for truncation of movements at the edge of the detector array.

If \code{RPSV} cannot be computed the algorithm tries to use observed
mean recapture distance \eqn{\bar{d}}{d-bar}.  Computation of
\eqn{\bar{d}}{d-bar} fails if there no recaptures, and all returned
values are NA.

If the mask has more than 100 points then a proportion 1--\code{thin} of
points are discarded at random to speed execution.

The argument \code{tol} is passed to \code{\link{uniroot}}.  It may be a
vector of two values, the first for g0 and the second for sigma.

If \code{traps(capthist)} has a \link{usage} attribute (defining effort
on each occasion at each detector) then the value of g0 is divided by
the mean of the non-zero elements of usage.  This adjustment is not
precise.

If \code{adjustg0} is TRUE then an adjustment is made to g0 depending
on the value of \code{binomN}.  For Poisson counts (\code{binomN = 0})
the adjustment is linear on effort (adjusted.g0 = g0 /
usage).  Otherwise, the adjustment is on the hazard scale (adjusted.g0 =
1 -- (1 -- g0) ^ (1 / (usage x binomN))).  An arithmetic average is taken
over all non-zero usage values (i.e. over used detectors and times). If
usage is not specified it is taken to be 1.0.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

}

\value{
  A list of parameter values :
  \item{ D }{Density (animals per hectare)}
  \item{ g0 }{Magnitude (intercept) of detection function }
  \item{ sigma }{Spatial scale of detection function (m) }
}

\note{
  
  \code{autoini} always uses the Euclidean distance between detectors and
  mask points.

  You may get this message from secr.fit: ``'autoini' failed to find g0;
  setting initial g0 = 0.1''.  If the fitted model looks OK (reasonable
  estimates, non-missing SE) there is no reason to worry about the
  starting values.  If you get this message and model fitting fails then
  supply your own values in the start argument of secr.fit.
  
}

\references{

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software for analysing capture--recapture data from passive detector arrays. \emph{Animal Biodiversity and Conservation} \bold{27},
 217--228.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\seealso{ \code{\link{capthist}}, \code{\link{mask}}, \code{\link{secr.fit}}, \code{\link{dbar}}  }

\examples{

\dontrun{

demotraps <- make.grid()
demomask <- make.mask(demotraps)
demoCH <- sim.capthist (demotraps, popn = list(D = 5, buffer = 100), seed = 321)
autoini (demoCH, demomask)

}

}
\keyword{ models }
