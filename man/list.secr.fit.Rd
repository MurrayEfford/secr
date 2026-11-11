\name{list.secr.fit}
\alias{list.secr.fit}

\title{Fit Multiple SECR Models}

\description{
  
  This function is a wrapper for \code{\link{secr.fit}} that allows multiple models to be fitted.

}

\usage{

list.secr.fit (..., constant = list(), prefix = "fit", names = NULL, ncores = NULL, bymodel = FALSE) 

}

\arguments{
  
  \item{\dots}{varying arguments of \code{secr.fit}}
  
  \item{constant}{list of named arguments held constant}
  
  \item{prefix}{character prefix for automatic names}
  
  \item{names}{character names of output}
  
  \item{ncores}{integer number of threads for secr.fit (bymodel = FALSE), or\cr 
      cores to be divided among worker processes (bymodel = TRUE)}
  
  \item{bymodel}{logical; if TRUE a separate worker process is used for each model}
  
}

\details{

The \dots argument may be one or several vectors of the same length that refer to a different named argument of \code{\link{secr.fit}}. \code{secr.fit} is called with the constant arguments plus the first value in each vector, then the second value, etc. The logic follows \code{\link{mapply}}.

Each of the \dots arguments may also be a named argument with a single value, although the compound values should be wrapped in list(), passed by name (in quotes), or placed in the 'constant' list to avoid misinterpretation. For example, the capthist argument of \code{secr.fit} should be be wrapped in list() or " " if it is placed outside 'constant'.

'prefix' is used only if 'names' is not supplied.

The argument `ncores' overrides ncores in ... or `constant' when `bymodel = FALSE'.

`bymodel = TRUE' with `ncores' > 1 causes a \pkg{\link{parallel}} cluster to be 
generated (PSOCK on Windows, FORK on Unix) and each model is fitted on a separate worker. 
If `ncores' exceeds the number of models specified then the cores are divided 
equally among the workers. 

If `ncores = NULL' then ncores takes the value from setNumThreads().
}

\value{

An \code{\link{secrlist}} of the successful model fits (see \code{\link{secr.fit}}). 

In the special case that `constant' specifies `details' with component LLonly = TRUE, 
the list of values from \code{secr.fit} without checking or modification.

From version 5.5.2 the total elapsed time in seconds is returned as attribute 'proctime'.
}

\seealso{
  \code{\link{secr.fit}},
  \code{\link{AIC.secr}},
  \code{\link{predict.secr}}
}

\note{

This function replaces the previous function \code{par.secr.fit}: since the 
introduction of multi-threading in \pkg{secr} 4.0 it is often no longer efficient 
to use parallel worker processes. However, the option of fitting each model on a 
separate worker process was added to \code{list.secr.fit} in version 5.5.2, 
restoring functionality from the deprecated function \code{par.secr.fit}.

In this version the random number seed is not set on the parallel workers, 
which may be a problem in rare cases when bymodel is used with an estimator 
that uses random numbers (e.g. nsim for simulating overdispersion).

}

\examples{

\dontrun{

# fit two detection models
fits <- list.secr.fit (model = c(g0~1, g0~b), constant = list(captdata, trace = FALSE))
AIC(fits)

# alternatively,
fits <- list.secr.fit ('captdata', model = c(g0~1, g0~b), trace = FALSE)
attr(fits,'proctime')
AIC(fits)

# using parallel cluster
fits <- list.secr.fit ('captdata', model = c(g0~1, g0~b), trace = FALSE, ncores=2, bymodel=TRUE)
attr(fits,'proctime')  # much slower
AIC(fits)

# replacing par.derived and par.region.N:

lapply(fits, derived)
lapply(fits, region.N)

}

}

