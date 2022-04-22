\name{makeStart}
\alias{makeStart}

\title{

Initial Parameter Values

}

\description{

\code{makeStart()} wraps the code previously internal to \code{secr.fit()}. 
It will not usually be called directly.

}

\usage{

makeStart(start = NULL, parindx, capthist, mask, detectfn, link, 
    details = NULL, fixed = NULL, CL = FALSE, anypoly = FALSE, 
    anytrans = FALSE, alltelem = FALSE, sighting = FALSE)

}

\arguments{
  \item{start}{optional starting values as described for \code{\link{secr.fit}}}
  \item{parindx}{list with one component per real parameter giving corresponding indices of coefficients (beta parameters) }
  \item{capthist}{ \code{\link{capthist}} object including capture data and detector (trap) layout }
  \item{mask}{ \code{\link{mask}} object or (for a multi-session analysis) a list of \code{mask} objects, one for each session}
  \item{detectfn}{integer code or character string for shape of detection function 0 = halfnormal, 1 = hazard rate etc. -- see
      \link{detectfn}}
  \item{link}{list of link function names ("log", "logit", "identity") for each real parameter }
  \item{details}{list of additional arguments (see \code{\link{secr.fit}})}
  \item{fixed}{list with optional components corresponding to real parameters giving the 
       scalar value to which the parameter is to be fixed }
  \item{CL}{logical, if true then model does not include density }
  \item{anypoly}{logical, TRUE if any polygon detectors}
  \item{anytrans}{logical, TRUE if any transect detectors}
  \item{alltelem}{logical, TRUE if any telemetry detectors}
  \item{sighting}{logical, TRUE if any sighting data}
  
}

\details{

If `start' is a previously fitted model then only `start' and `parindx' are required.

}

\value{

Numeric vector with one value for each coefficient (beta parameter) in model.

}

\seealso{

\code{\link{secr.fit}}

}

\examples{

makeStart(secrdemo.0, list(D = 1, g0 = 2:3, sigma = 4))

}

