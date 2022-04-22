\name{timevaryingcov}
\alias{timevaryingcov}
\alias{timevaryingcov<-}

\title{ Time-varying Covariates }

\description{
  Extract or replace time varying covariates
}

\usage{
timevaryingcov(object, ...)
timevaryingcov(object) <- value
}

\arguments{
  \item{object}{ an object of class \code{traps} or \code{capthist}}
  \item{value}{ a list of named vectors }
  \item{\dots}{ other arguments (not used) }
}

\details{

  The timevaryingcov attribute is a list of one or more named vectors. 
  Each vector identifies a subset of columns of covariates(object), one for each occasion. If character values are used they should correspond to covariate names.

In \pkg{secr} models, time-varying covariates are restricted to traps objects.  
Time-varying (session-specific) individual covariates may be used in \pkg{openCR}. 
The following remarks apply to time-varying traps covariates.

  The name of the vector may be used in a model formula; when the model
  is fitted, the value of the trap covariate on a particular occasion is
  retrieved from the column indexed by the vector.
  
  For replacement, if \code{object} already has a \code{\link{usage}}
  attribute, the length of each vector in \code{value} must match exactly
  the number of columns in \code{usage(object)}.
  
  When converting a multi-session capthist object into a robust-design ``single-session'' object with function \code{join} the argument `timevaryingcov' is used to collate 
  covariate values across sessions in a form suitable for inclusion in 
  \pkg{openCR} models (see \code{\link{join}}).
  
}

\value{
  \code{timevaryingcov(object)} returns the timevaryingcov attribute of
  \code{object} (may be NULL).
}

\seealso{

\code{\link{join}}

}

\note{

It is usually better to model varying effort directly, via the \link{usage}
attribute (see \href{https://www.otago.ac.nz/density/pdfs/secr-varyingeffort.pdf}{secr-varyingeffort.pdf}).

Models for data from detectors of type `multi', `polygonX' or `transectX' take much
longer to fit when detector covariates of any sort are used.

Time-varying covariates are not available with the (default) 'fastproximity' option.

See \href{https://www.otago.ac.nz/density/pdfs/secr-varyingeffort.pdf}{secr-varyingeffort.pdf} for input of detector covariates from a file.

}

\examples{

# make a trapping grid with simple covariates
temptrap <- make.grid(nx = 6, ny = 8, detector = "multi") 
covariates (temptrap) <- data.frame(matrix(
    c(rep(1,48*3),rep(2,48*2)), ncol = 5))
head(covariates (temptrap))

# identify columns 1-5 as daily covariates
timevaryingcov(temptrap) <- list(blockt = 1:5)
timevaryingcov(temptrap)

\dontrun{

# default density = 5/ha, noccasions = 5
CH <- sim.capthist(temptrap, detectpar = list(g0 = c(0.15, 0.15,
    0.15, 0.3, 0.3), sigma = 25))

fit.1 <- secr.fit(CH, trace = FALSE) 
fit.tvc2 <- secr.fit(CH, model = g0 ~ blockt, trace = FALSE) 

# because variation aligns with occasions, we get the same with:
fit.t2 <- secr.fit(CH, model = g0 ~ tcov, timecov = c(1,1,1,2,2),
    trace = FALSE) 

predict(fit.t2, newdata = data.frame(tcov = 1:2))
predict(fit.tvc2, newdata = data.frame(blockt = 1:2))

# now model some more messy variation
covariates (traps(CH))[1:10,] <- 3
fit.tvc3 <- secr.fit(CH, model = g0 ~ blockt, trace = FALSE) 

AIC(fit.tvc2, fit.t2, fit.tvc3)
# fit.tvc3 is the 'wrong' model

}

}

\keyword{ manip }

