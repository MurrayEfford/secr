\name{list.secr.fit}
\alias{list.secr.fit}

\title{Fit Multiple SECR Models}

\description{
  
  This function is a wrapper for \code{\link{secr.fit}} that allows multiple models to be fitted.

}

\usage{

list.secr.fit (..., constant = list(), prefix = "fit", names = NULL) 

}

\arguments{
  
  \item{\dots}{varying arguments of \code{secr.fit}}
  
  \item{constant}{list of named arguments held constant}
  
  \item{prefix}{character prefix for automatic names}
  
  \item{names}{character names of output}
  
}

\details{

The \dots argument may be one or several vectors of the same length that refers to a different named argument of \code{\link{secr.fit}}. \code{secr.fit} is called with the constant arguments plus the first value in each vector, then the second value, etc. The logic follows \code{\link{mapply}}.

Each of the \dots arguments may also be a named argument with a single value, although the compound values should be wrapped in list(), passed by name (in quotes), or placed in the 'constant' list to avoid misinterpretation. For example, the capthist argument of \code{secr.fit} should be be wrapped in list() or " " if it is placed outside 'constant'.

'prefix' is used only if 'names' is not supplied.

}

\value{

An \code{\link{secrlist}} of model fits (see \code{\link{secr.fit}}). 

}

\seealso{
  \code{\link{secr.fit}},
  \code{\link{AIC.secr}},
  \code{\link{predict.secr}}
}

\note{

This function replaces the previous function \code{par.secr.fit}: since the introduction of multi-threading in \pkg{secr} 4.0 it is no longer efficient to use parallel worker processes.

}

\examples{

\dontrun{

# fit two detection models
fits <- list.secr.fit (model = c(g0~1, g0~b), constant = list(captdata, trace = FALSE))
AIC(fits)

# alternatively,
fits <- list.secr.fit ('captdata', model = c(g0~1, g0~b), trace = FALSE)
AIC(fits)

# replacing par.derived and par.region.N:

lapply(fits, derived)
lapply(fits, region.N)

}

}

