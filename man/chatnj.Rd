\name{chat}
\alias{chat.nj}
\alias{adjustVarD}

\title{ Overdispersion of Activity Centres}

\description{

Stochastic variation in the intensity surface may cause the number of detected 
individuals \eqn{n} to be overdispersed relative to a Poisson distribution 
(Efford and Fletcher 2025). This can cause the sampling variance of density 
estimates to be understated. 

Use \code{chat.nj} to compute Fletcher's \eqn{\hat c} estimate of overdispersion 
for use as a variance inflation factor. The method requires replicate samples of 
the intensity surface, such as from multiple independent subgrids. 
Replicates may be provided as sessions or clusters of detectors within a session 
(see \link{cluster} and \code{\link{trap.builder}}).

\code{adjustVarD} adjusts the SE and confidence limits of density estimates 
using the given \eqn{\hat c}. The implementation is limited to simple detection models 
(see Warnings).

See Cooch and White (2022) for an introduction to measurement of overdispersion 
in capture--recapture due to non-independence in the detection process. 
The focus here is on overdispersion of the number detected \eqn{n} relative to 
the Poisson (or binomial) distribution of \eqn{n} expected when the 
distribution of activity centres is an inhomogeneous or homogeneous Poisson (or binomial) 
point process (IHPP).

}

\usage{

chat.nj(object, bycluster = FALSE, ...)

adjustVarD(object, chatmin = 1, alpha = 0.05, chat = NULL)

}

\arguments{

  \item{object}{fitted secr model or dataframe (see Warnings for restrictions)}
  \item{bycluster}{logical; if TRUE then j refers to clusters within a 
      single-session capthist rather than sessions}
  \item{\dots}{other arguments passed to \code{\link{Fletcher.chat}}, e.g., verbose = FALSE}
  
  \item{chatmin}{minimum value of Fletcher's \eqn{\hat c}}
  \item{alpha}{alpha level for confidence intervals }
  \item{chat}{numeric chat (optional)}
}

\details{

\code{chat.nj} uses \code{\link{expected.n}} to compute the expected values.

No adjustment is made by \code{adjustVarD} when \eqn{\hat c} is less than the 
minimum. Set 'chatmin' to zero to override.

\code{adjustVarD} also accepts a single dataframe as the argument `object';
the dataframe should have row `D' and columns `link', `estimate', `SE.estimate' 
as in the output from \code{\link{predict.secr}}.

}

\value{

For \code{chat.nj}, usually a list comprising --

\item{expected}{expected number at each cluster or session} 
\item{observed}{observed number at each cluster or session}
\item{stats}{vector of summary statistics: mean(expected), var(expected), 
mean(observed), var(observed), sbar, nu (=df), cX2 = X2/nu}
\item{chat}{\eqn{\hat c} (Fletcher or Wedderburn depending on `type')}

If `verbose = FALSE' then only the numeric value of \eqn{\hat c} is returned (a vector of 2 values if `type = "both"').

For \code{adjustVarD}, a dataframe with one row for each session, based on \code{\link{predict.secr}} or \code{\link{derived.secr}}, with extra column `c-hat'.

}

\section{Warning}{

The variance inflation factor given by \code{chat.nk} was shown by Efford and 
Fletcher (2025) to be inadequate and should not be used. For replicate spatial 
samples, \code{\link{chat.nj}} is a better alternative.

\code{adjustVarD} previously computed Fletcher's `chat' using \code{chat.nk}; 
that is no longer recommended.

}

\seealso{ 

  \code{\link{secr}}, 
  \code{\link{make.mask}},
  \code{\link{Detection functions}}, 
  \code{\link{Fletcher.chat}},
  \code{\link{cluster}}

}
  
\references{

Bischof, R., P. Dupont, C. Milleret, J. Chipperfield, and J. A. Royle. 2020. 
Consequences of ignoring group association in spatial capture--recapture analysis. 
\emph{Wildlife Biology} wlb.00649. \doi{10.2981/wlb.00649}

Cooch, E. and White, G. (eds) (2022) \emph{Program MARK: A Gentle Introduction}. 
  22nd edition. Most recent edition available online at www.phidot.org/software/mark/docs/book/.

Efford, M. G. and D. Fletcher. 2025. Effect of spatial overdispersion on 
confidence intervals for population density estimated by spatial capture-recapture.
\emph{bioRxiv} https://doi.org/10.1101/2024.03.12.584742 

Fletcher, D. (2012) Estimating overdispersion when fitting a generalized linear model to sparse data.
\emph{Biometrika} \bold{99}, 230--237.

Wedderburn, R. W. M. (1974) Quasi-likelihood functions, generalized linear models, and the Gauss-Newton
method. \emph{Biometrika} \bold{61}, 439--47.

}

\examples{

\dontrun{

## Clustered design
mini <- make.grid(nx = 3, ny = 3, spacing = 50, detector = "proximity")

tempgrids2 <- trap.builder (
    cluster = mini , method = "all",
    frame = expand.grid(x = seq(200, 1800, 400), y = seq(200, 1800, 400)), 
    plt = TRUE)
tempmask2 <- make.mask(tempgrids2, buffer = 100, spacing = 20)
pop1 <- sim.popn(D = 2, core = tempmask2, model2D = 'rLGCP', 
                 details = list(var = 0.5, scale = 100))
capt <- sim.capthist(tempgrids2, popn = pop1, noccasions = 5,
    detectpar = list(g0 = 0.2, sigma = 25))

fit <- secr.fit(capt, mask = tempmask2, trace = FALSE)

chat.nj(fit, bycluster = TRUE)

}

}


