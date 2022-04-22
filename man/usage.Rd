\name{usage}
\alias{usage}
\alias{usage<-}
\alias{effort}

\title{ Detector Usage }

\description{
  Extract or replace usage (effort) information of a \code{traps} object
  (optional).
}

\usage{
usage(object, \dots)
usage(object) <- value
}
\arguments{
  \item{object}{ \code{traps} object }
  \item{value}{ numeric matrix of detectors x occasions }
  \item{\dots}{ other arguments (not used) }
}
\details{
  
  In secr versions before 2.5.0, usage was defined as a binary value
  (1 if trap \eqn{k} used on occasion \eqn{s}, zero otherwise).

  In later versions, usage may take nonnegative real values and will be
  interpreted as effort. This corresponds to the constant T_s used for
  the duration of sampling by Borchers and Efford (2008). Effort is
  modelled as a known linear coefficient of detection probability on the
  hazard scale
  (\href{https://www.otago.ac.nz/density/pdfs/secr-varyingeffort.pdf}{secr-varyingeffort.pdf};
  Efford et al. 2013).
  
  For replacement of usage, various forms are possible for \code{value}:

  - a matrix in which the number of rows of \code{value} exactly
  matches the number of traps \emph{K} in \code{object}

  - a vector of two values, the usage (typically 1) and the number of
  occasions \emph{S} (a \emph{K} x \emph{S} matrix will be filled with
  the first value)

  - a vector of \emph{R+1} values where \emph{R} is the number of
  sessions in a multi-session object and elements \emph{2..R+1} correspond to
  the numbers of occasions \emph{S1, S2,...} in each session

  - the usage only (typically 1) (only works when replacing an existing
  usage matrix with known number of occasions).  
  
}
\value{
  
  usage(object) returns the usage matrix of the \code{traps}
  object. \code{usage(object)} may be NULL.
  
}

\references{
Efford, M. G., Borchers D. L. and Mowat, G. (2013) Varying effort in
  capture--recapture studies. \emph{Methods in Ecology
    and Evolution} \bold{4}, 629--636.
}

\note{

  At present, assignment of usage to the traps objects of a multisession
  capthist object results in the loss of session names from the latter.
  
}

\seealso{
  
  \code{\link{traps}},
  \code{\link{usagePlot}},
  \code{\link{read.capthist}},
  \code{\link{addSightings}}
}

\examples{
demo.traps <- make.grid(nx = 6, ny = 8)
## random usage over 5 occasions
usage(demo.traps) <- matrix (sample(0:1, 48*5, replace = TRUE, 
    p = c(0.5,0.5)), nc = 5)
usage(demo.traps)
summary(demo.traps)

usage(traps(ovenCH)) <- c(1,9,10,10,10,10)
## restore lost names
names(ovenCH) <- 2005:2009

}
\keyword{ manip }

