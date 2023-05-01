\name{chat}
\alias{Enk}
\alias{chat.nk}
\alias{adjustVarD}

\title{ Overdispersion of Activity Centres}

\description{

Activity centres may be clumped (overdispersed) relative to a Poisson distribution,
the model used in \code{\link{secr.fit}} (Borchers and Efford 2008). This causes 
the sampling variance of density estimates to be understated. The problem is 
addressed with a variance inflation factor, essentially a measure of overdispersion.

Functions described here compute the expected number of individuals detected at each detector 
(Enk) and use that to compute Fletcher's \eqn{\hat c} estimate of overdispersion 
for use as a variance inflation factor. 

\code{Enk} uses exact formulae for 'multi', 'proximity' and 'count' detector types. Other types may be simulated by setting a positive value for 'nrepl', which should be large (e.g., nrepl = 10000).

\code{adjustVarD} adjusts the SE and confidence limits of density estimates 
using Fletcher's \eqn{\hat c}. The implementation is limited to simple detection models 
(see Warnings).

See Cooch and White (2022) for an introduction to measurement of overdispersion 
in capture--recapture. The focus here is on overdispersion of activity centres 
relative to a Poisson distribution, rather than on non-independence in the 
spatial detection process.

}

\usage{

Enk(D, mask, traps, detectfn = 0, detectpar = list(g0 = 0.2,
    sigma = 25, z = 1), noccasions = NULL, binomN = NULL,
    userdist = NULL, ncores = NULL, nrepl = NULL) 

chat.nk(object, nsim = NULL, ...)

adjustVarD(object, chatmin = 1, alpha = 0.05, chat = NULL)

}

\arguments{
  \item{D}{numeric density, either scalar or vector of length nrow(mask)}
  \item{mask}{single-session habitat mask}
  \item{traps}{ \code{traps} object }
  \item{detectfn}{ integer code for detection function q.v. }
  \item{detectpar}{ a named list giving a value for each parameter of detection function }
  \item{noccasions}{ number of sampling intervals (occasions) }
  \item{binomN}{ integer code for discrete distribution (see
    \code{\link{secr.fit}}) }
  \item{userdist}{user-defined distance function or matrix (see \link{userdist})}
  \item{ncores}{integer number of threads}
  \item{nrepl}{integer number of replicates for E(nk) by simulation (optional)}

  \item{object}{fitted secr model or dataframe (see Warnings for restrictions)}
  \item{nsim}{integer number of c-hat values to simulate (optional) }
  \item{...}{other arguments passed to \code{\link{Fletcher.chat}} (verbose, type)}

\item{chatmin}{minimum value of Fletcher's \eqn{\hat c}}

\item{alpha}{alpha level for confidence intervals }
\item{chat}{numeric chat (optional)}
}

\details{

If \code{traps} has a \link{usage} attribute then \code{noccasions} is
set accordingly; otherwise it must be provided.

The environment variable RCPP_PARALLEL_NUM_THREADS determines the number of 
parallel threads. It is set to the value of \code{ncores}, unless that is NULL 
(see \code{\link{setNumThreads}}).

The conventional lack-of-fit statistic is
\eqn{\hat c_X = X^2/(K-p)} where \eqn{K} is the number of detectors,  \eqn{p} is the number of estimated parameters, and 
\deqn{X^2 = \sum_k (n_k - E (n_k))^2/  E(n_k).} 

Fletcher's \eqn{\hat c} is an improvement on \eqn{\hat c_X} that is less affected 
by small expected counts. It is defined by
\deqn{\hat c = c_X / (1+ \bar s),} 
where \eqn{\bar s = \sum_k s_k / K} and \eqn{s_k = (n_k - E(n_k)) / E(n_k)}.
 
No adjustment is made by \code{adjustVarD} when \eqn{\hat c} is less than the 
minimum. \code{adjustVarD} by default computes Fletcher's `chat' using \code{chat.nk}, 
but a value may be provided.

If chat has been computed separately and provided in the argument of that name, \code{adjustVarD} also accepts a single dataframe as the argument `object'; the dataframe should have row `D' and columns `link', `estimate', `SE.estimate' as in the output from \code{\link{predict.secr}}.

}

\value{

For \code{Enk}, a vector of expected counts, one for each detector in \code{traps}.

For \code{chat.nk}, a list comprising --

\item{expected.nk}{expected number at each detector} 
\item{nk}{observed number at each detector}
\item{stats}{vector of summary statistics: mean(expected.nk), var(expected.nk), mean(nk), var(nk), nu (=df), X2/nu}
\item{chat}{Fletcher's \eqn{\hat c}}
 
For \code{adjustVarD}, a dataframe with one row for each session, based on \code{\link{predict.secr}} or \code{\link{derived.secr}}, with extra column `c-hat'.

}

\section{Warning}{

These functions are experimental in secr 4.5.11, and do not work with 
polygon-like and single-catch detectors. No allowance is made for modelled variation in 
detection parameters with respect to occasion, detector or animal; 
this includes mixture models (e.g., g0~h2).

Versions before 4.5.11 did not correctly compute expected counts for multi-catch detectors. 
}

\seealso{ 

  \code{\link{secr}}, 
  \code{\link{make.mask}},
  \code{\link{Detection functions}}, 
  \code{\link{Fletcher.chat}}

}
  
\references{

Cooch, E. and White, G. (eds) (2022) \emph{Program MARK: A Gentle Introduction}. 
  22nd edition. Available online at http://www.phidot.org/software/mark/docs/book/.

Fletcher, D. (2012) Estimating overdispersion when fitting a generalized linear model to sparse data.
\emph{Biometrika} \bold{99}, 230--237.

}

\examples{
\donttest{
  temptrap <- make.grid()
  msk <- make.mask(temptrap)
  ## expected number of individuals per detector (multi-catch) 
  Enk (D = 5, msk, temptrap, detectpar = list(g0 = 0.2, sigma = 25),
    noccasions = 5)

}

}

\keyword{ manip }

