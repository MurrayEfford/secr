\name{write.captures}
\alias{write.captures}
\alias{write.traps}
\alias{write.mask}

\title{ Write Data to Text File }

\description{
  Export detections or detector layout or mask to a text file in format suitable
  for input to DENSITY. 
}

\usage{

write.captures(object, file = "", deblank = TRUE, header = TRUE,
    append = FALSE, sess = "1", ndec = 2, covariates = FALSE, tonumeric
    = TRUE, ...)

write.traps(object, file = "", deblank = TRUE,  header = TRUE,
    ndec = 2, covariates = FALSE, ...)

write.mask(object, file = "", header = TRUE, ndec = 2, covariates = TRUE, ...)

}
\arguments{
  \item{object}{ \code{capthist} or \code{traps} object }
  \item{file}{ character name of output file }
  \item{deblank}{ logical; if TRUE remove any blanks from character string used to identify detectors }
  \item{header}{ logical; if TRUE output descriptive header }
  \item{append}{ logical; if TRUE output is appended to an existing file }
  \item{sess}{ character session identifier }
  \item{ndec}{ number of digits after decimal point for x,y coordinates }
  \item{covariates}{ logical or a character vector of covariates to
    export }
  \item{tonumeric}{ logical for whether factor and character covariates should be
    converted to numeric values on output}
  \item{\dots}{ other arguments passed to \code{write.table} }
}
\details{
  
  Existing file will be replaced without warning if \code{append =
  FALSE}. In the case of a multi-session capthist file, session names
  are taken from \code{object} rather than \code{sess}.

  \code{\link{write.capthist}} is generally simpler to use if you want to export
  both the capture data and trap layout from a \code{capthist} object.

  By default individual covariates are not exported. When exported they
  are repeated for each detection of an individual. Factor covariates
  are coerced to numeric before export.
  
  For \code{write.mask}, \code{header = TRUE} also causes column names to be exposed.
  
}
\value{

None

}
\seealso{
\code{\link{as.data.frame.capthist}}
}

\examples{

  write.captures (captdata)

}
\keyword{ IO }
