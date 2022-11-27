\name{chat}
\alias{Enk}
\alias{chat.nk}

\title{ Overdispersion of Activity Centres}

\description{

Activity centres may be clumped (overdispersed) relative to a Poisson distribution,
the model used in \code{\link{secr.fit}} (Borchers and Efford 2008). This causes 
the sampling variance of density estimates to be understated. The problem is 
addressed with a variance inflation factor, essentially a measure of overdispersion.

Functions described here compute the expected number of individuals detected at each detector 
(Enk) and use that to compute Fletcher's \eqn{\hat c} estimate of overdispersion 
for use as a variance inflation factor. 

See Cooch and White (2022) for an introduction to measurement of overdispersion 
in capture--recapture. The focus here is on overdispersion of activity centres 
relative to a Poisson distribution, rather than on non-independence in the 
spatial detection process.

}

\usage{

Enk(D, mask, traps, detectfn = 0, detectpar = list(g0 = 0.2,
    sigma = 25, z = 1), noccasions = NULL, binomN = NULL,
    userdist = NULL, ncores = NULL) 

chat.nk(object)

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

\item{object}{fitted secr model (see Warnings for restrictions)}

}
\details{

If \code{traps} has a \link{usage} attribute then \code{noccasions} is
set accordingly; otherwise it must be provided.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

The conventional lack-of-fit statistic is
\eqn{\hat c_X = X^2/(K-1)} where \eqn{K} is the number of detectors and 
\deqn{X^2 = \sum_k (n_k - E (n_k))^2/  E(n_k).} 

Fletcher's \eqn{\hat c} is an improvement on \eqn{\hat c_X} that is less affected by sparse counts. 
It is defined by
\deqn{\hat c = c_X / (1+ \bar s),} 
where \eqn{\bar s = \sum_k s_k / K} and \eqn{s_k = (n_k - E(n_k)) / E(n_k)}.
 
}

\value{

For \code{Enk}, a vector of expected counts, one for each detector in \code{traps}.

For \code{chat.nk}, a list comprising --

\item{expected.nk}{expected number at each detector} 
\item{nk}{observed number at each detector}
\item{stats}{vector of summary statistics: mean(expected.nk), var(expected.nk), mean(nk), var(nk), X2/nu}
\item{chat}{Fletcher's \eqn{\hat c}}
 
}

\section{Warning}{

These functions are experimental in secr 4.5.8, and do not work with 
polygon-like detectors. No allowance is made for varying detection parameters
with respect to occasion, detector, animal.

}

\seealso{ 

  \code{\link{secr}}, 
  \code{\link{make.mask}},
  \code{\link{Detection functions}}, 
  \code{\link{pdot}}

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
  ## expected number of individuals per detector. 
  Enk (D = 5, msk, temptrap, detectpar = list(g0 = 0.2, sigma = 25),
    noccasions = 5)

  chat.nk(secrdemo.0)    
}

}

\keyword{ manip }

