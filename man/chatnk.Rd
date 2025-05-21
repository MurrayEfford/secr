\name{nk}
\alias{nk}
\alias{Enk}
\alias{chat.nk}

\title{ Number per detector}

\description{

Functions described here compute the observed (nk) or expected (Enk) number of 
individuals detected at each detector. 

The function chat.nk computes an erroneous estimate of the overdispersion of 
activity centres (Bischof et al. 2020; Efford and Fletcher 2025).

\code{Enk} uses exact formulae for 'multi', 'proximity' and 'count' detector types. 
Other types may be simulated by setting a positive value for 'nrepl', which 
should be large (e.g., nrepl = 10000).


}

\usage{

nk(capthist)

Enk(D, mask, traps, detectfn = 0, detectpar = list(g0 = 0.2,
    sigma = 25, z = 1), noccasions = NULL, binomN = NULL,
    userdist = NULL, ncores = NULL, nrepl = NULL) 

chat.nk(object, nsim = NULL, ...)

}

\arguments{
  \item{capthist}{secr \code{capthist} object}
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

}

\details{

If \code{traps} has a \link{usage} attribute then \code{noccasions} is
set accordingly; otherwise it must be provided.

The environment variable RCPP_PARALLEL_NUM_THREADS determines the number of 
parallel threads. It is set to the value of \code{ncores}, unless that is NULL 
(see \code{\link{setNumThreads}}).

\code{chat.nk} may be used to simulate \eqn{\hat c} values under the given 
model (set nsim > 0). The \dots argument may include 'ncores = x' (x>1) to 
specify parallel processing of simulations - the speed up is large on unix-like 
machines for which the cluster type of \code{\link{makeCluster}} is "FORK" 
rather than "PSOCK". If 'ncores' is not provided then the value returned by 
\code{setNumThreads()} is used.

}

\value{

For \code{nk}, a vector of observed counts, one for each detector in \code{traps(capthist)}.

For \code{Enk}, a vector of expected counts, one for each detector in \code{traps}.

For \code{chat.nk}, usually a list comprising --

\item{expected.nk}{expected number at each detector} 
\item{nk}{observed number at each detector}
\item{stats}{vector of summary statistics: mean(expected.nk), var(expected.nk), mean(nk), var(nk), nu (=df), X2/nu}
\item{chat}{\eqn{\hat c} (Fletcher or Wedderburn depending on `type')}
There are two variations --
 
If `verbose = FALSE' then only the numeric value of \eqn{\hat c} is returned (a vector of 2 values if `type = "both"').

If \code{chat.nk} is called with `nsim > 0' then the output is a list comprising --

\item{type}{from input}
\item{nsim}{from input}
\item{sim.chat}{vector of simulated \eqn{\hat c}} 
\item{chat}{\eqn{\hat c} (Fletcher or Wedderburn depending on `type')}
\item{p}{probability of observing \eqn{\hat c} as large as this (from rank of chat among sim.chat)}

}

\section{Warning}{

The variance inflation factor given by \code{chat.nk} was shown by Efford and 
Fletcher (2025) to be inadequate and should not be used. For replicate spatial 
samples, \code{\link{chat.nj}} is a better alternative.

These functions are experimental and do not work with 
polygon-like and single-catch detectors. No allowance is made for modelled variation in 
detection parameters with respect to occasion, detector or animal; 
this includes mixture models (e.g., g0~h2).

Versions before 4.5.11 did not correctly compute expected counts for multi-catch detectors. 

}

\seealso{ 

  \code{\link{secr}}, 
  \code{\link{make.mask}},
  \code{\link{Detection functions}}, 
  \code{\link{Fletcher.chat}},
  \code{\link{chat.nj}}

}
  
\references{

Bischof, R., P. Dupont, C. Milleret, J. Chipperfield, and J. A. Royle. 2020. 
Consequences of ignoring group association in spatial capture--recapture analysis. 
\emph{Wildlife Biology} wlb.00649. \doi{10.2981/wlb.00649}

Efford, M. G. and D. Fletcher. 2025. Effect of spatial overdispersion on 
confidence intervals for population density estimated by spatial capture-recapture.
\emph{bioRxiv} https://doi.org/10.1101/2024.03.12.584742 

Fletcher, D. (2012) Estimating overdispersion when fitting a generalized linear model to sparse data.
\emph{Biometrika} \bold{99}, 230--237.

Wedderburn, R. W. M. (1974) Quasi-likelihood functions, generalized linear models, and the Gauss-Newton
method. \emph{Biometrika} \bold{61}, 439--47.

}

\examples{
\donttest{
  temptrap <- make.grid()
  msk <- make.mask(temptrap)
  ## expected number of individuals per detector (multi-catch) 
  Enk (D = 5, msk, temptrap, detectpar = list(g0 = 0.2, sigma = 25),
    noccasions = 5)

}

# useful plotting function for simulated chat (nsim>0)
plotchat <- function(chat, head = '', breaks = seq(0.5,2.5,0.05)) {
    hist(chat$sim.chat, xlim = range(breaks), main = head, xlab = 'c-hat',
        breaks = breaks, cex.main = 1, yaxs = 'i')
    abline(v = chat$chat, lwd = 1.5, col = 'blue')
}

}


