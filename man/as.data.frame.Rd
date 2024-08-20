\name{as.data.frame}
\alias{as.data.frame}
\alias{as.array}
\alias{as.data.frame.capthist}
\alias{as.array.capthist}
\alias{as.data.frame.traps}

\title{ Coerce capthist to Data Frame or Array }

\description{

  Method for generic \code{as.data.frame} function that partially reverses \code{make.capthist}.

}

\usage{

\method{as.data.frame}{capthist}(x, row.names = NULL, optional = FALSE, covariates = FALSE, 
                                    fmt = c("trapID", "XY"), ...)
\method{as.data.frame}{traps}(x, row.names = NULL, optional = FALSE, usage = FALSE, 
                              covariates = FALSE, ...)
\method{as.array}{capthist}(x, ...)

}

\arguments{
  \item{x}{ \code{capthist}  object }
  \item{row.names}{unused argument of generic function}
  \item{optional}{unused argument of generic function}
  \item{covariates}{ logical or a character vector of covariates to export }
  \item{fmt}{ character string for capture format }
  \item{usage}{ logical; if TRUE then usage columns are appended if present }
  \item{\dots}{ other arguments (not used) }
  
}

\details{
  
  By default individual covariates are not exported. When exported they
  are repeated for each detection of an individual.
  
}
\value{

A data frame or list of data frames (in the case of a multisession input). 

For capthist objects --

The core columns are (Session, ID, Occasion, TrapID) or (Session, ID, Occasion, x, y), 
depending on the value of \code{fmt}.  Additional columns for covariates and signal 
strength (detector `signal') are appended to the right. 

For traps objects --

The core columns are (x, y). Usage columns are named u1, u2, ..., uS where S is the number of occasions.

The \code{as.array} method for capthist objects returns an object with the same dimensions and dimnames but different class, or a list of such objects in the case of multisession input.

}

\examples{

  as.data.frame (captdata)
  as.data.frame (traps(captdata))
  
}
\keyword{ IO }
