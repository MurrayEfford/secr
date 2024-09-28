\name{secr.design.MS}
\alias{secr.design.MS}
\alias{make.lookup}
\alias{insertdim}

\title{ Construct Detection Model Design Matrices and Lookups }

\description{
  Internal functions used by \code{\link{secr.fit}}. 
}

\usage{

secr.design.MS (capthist, models, timecov = NULL, sessioncov = NULL, 
    groups = NULL, hcov = NULL, dframe = NULL, naive = FALSE, CL = FALSE,
    keep.dframe = FALSE, full.dframe = FALSE, ignoreusage = FALSE, contrasts = NULL, ...)

make.lookup (tempmat)

insertdim (x, dimx, dims) 

}
\arguments{
  \item{capthist}{ \code{capthist} object }
  \item{models}{ list of formulae for parameters of detection }
  \item{timecov}{ optional dataframe of values of time (occasion-specific) covariate(s). }
  \item{sessioncov}{ optional dataframe of values of session-specific covariate(s). }
  \item{groups}{ optional vector of one or more variables with which to
    form groups. Each element should be the name of a factor variable in
    the \code{covariates} attribute of \code{capthist}. }
  \item{hcov}{ character name of an individual (capthist) covariate for known class
    membership in h2 models}
  \item{dframe}{ optional data frame of design data for detection parameters }
  \item{naive}{ logical if TRUE then modelled detection probability is for a naive animal (not caught previously); if FALSE then detection probability is contingent on individual's history of detection }
  \item{CL}{ logical; TRUE for model to be fitted by maximizing the conditional likelihood }
  \item{keep.dframe}{logical; if TRUE the dataframe of design data is
    included in the output}
  \item{full.dframe}{logical; if FALSE then padding rows are purged from
    output dframe (ignored if \code{keep.dframe = FALSE})}
  \item{ignoreusage}{logical; if TRUE any usage attribute of traps(capthist)
  is ignored}
  \item{contrasts}{contrast specification as for \code{\link[stats]{model.matrix}}}
  \item{\dots}{ other arguments passed to the \R function
    \code{\link{model.matrix}} }
  \item{tempmat}{matrix for which row lookup required}
  
  \item{x}{vector of character, numeric or factor values }
  \item{dimx}{vector of notional dimensions for x to fill in target array}
  \item{dims}{vector of notional dimensions of target array}
  
}
\details{

These are internal \pkg{secr} functions that you are unlikely ever to
use. 

Each real parameter is notionally different for each unique combination
of session, individual, occasion, detector and latent class, i.e., for \eqn{R}
sessions, \eqn{n} individuals, \eqn{S} occasions, \eqn{K} detectors and \eqn{M} latent classes 
there are \emph{potentially} \eqn{R \times n \times S \times K \times M}{R x n x
S x K x M} different values. Actual models always predict a \emph{much} reduced set
of distinct values, and the number of rows in the design matrix is
reduced correspondingly; a parameter index array allows these to
retrieved for any combination of session, individual, occasion and
detector.

The \code{keep.dframe} option is provided for the rare occasions that a
user may want to check the data frame that is an intermediate step in
computing each design matrix with \code{\link{model.matrix}} (i.e. the
data argument of \code{model.matrix}).

\dots may be used to pass \code{contrasts.arg} to \code{model.matrix}. 

}
\value{
  
  For \code{secr.design.MS}, a list with the components
  
  \item{designMatrices}{ list of reduced design matrices, one for each
    real detection parameter }
  
  \item{parameterTable}{ index to row of the reduced design matrix for
  each real detection parameter; dim(parameterTable) = c(uniquepar, np),
  where uniquepar is the number of unique combinations of paramater
  values (uniquepar < \eqn{RnSKM}) and np is the number of parameters in
  the detection model. }

  \item{PIA}{ Parameter Index Array - index to row of parameterTable for
    a given session, animal, occasion and detector; dim(PIA) =
    c(R,n,S,K,M) }

  \item{R}{ number of sessions}

  If \code{models} is empty then all components are NULL except for PIA
  which is an array of 1's (M set to 1).  

  Optionally (\code{keep.dframe = TRUE}) -
  \item{dframe}{ dataframe of design data, one column per covariate, one
    row for each c(R,n,S,K,M). For multi-session models n, S, and K
  refer to the maximum across sessions}
  \item{validdim}{ list giving the valid dimensions (n, S, K, M) before padding}

  For \code{make.lookup}, a list with components
  \item{lookup}{ matrix of unique rows}
  \item{index}{ indices in lookup of the original rows}
  
  For \code{insertdim}, a vector with length prod(dims) containing the 
  values replicated according to dimx.
  
}

\examples{
secr.design.MS (captdata, models = list(g0 = ~b))$designMatrices
secr.design.MS (captdata, models = list(g0 = ~b))$parameterTable

## peek at design data constructed for learned response model
head(captdata)
temp <- secr.design.MS (captdata, models = list(g0 = ~b),
    keep.dframe = TRUE)
a1 <- temp$dframe$animal == 1 &  temp$dframe$detector \%in\% 8:10
temp$dframe[a1,]

## ... and trap specific learned response model
temp <- secr.design.MS (captdata, models = list(g0 = ~bk),
    keep.dframe = TRUE)
a1 <- temp$dframe$animal == 1 &  temp$dframe$detector \%in\% 8:10
temp$dframe[a1,]

## place values 1:6 in different dimensions
insertdim(1:6, 1:2, c(2,3,6))
insertdim(1:6, 3, c(2,3,6))

}

\seealso{ \code{\link{D.designdata}}, \code{\link{model.matrix}} }

