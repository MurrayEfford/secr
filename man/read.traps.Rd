\name{read.traps}
\alias{read.traps}
\title{ Read Detector Data From File }

\description{
  Construct an object of class \code{traps} with detector locations from a text file or data frame. Usage per occasion and covariates may be included. Data may also be read from an Excel spreadsheet (see \href{../doc/secr-datainput.pdf}{secr-datainput.pdf}).

}

\usage{
read.traps(file = NULL, data = NULL, detector = "multi", covnames =
NULL, binary.usage = TRUE, markocc = NULL, trapID = NULL, ...)
}

\arguments{
  \item{file}{ character string with name of text file }
  \item{data}{ data frame of detector coordinates }
  \item{detector}{ character string for detector type }
  \item{covnames}{ character vector of names for detector covariate
    fields }
  \item{binary.usage}{logical; if FALSE will read usage fields as
  continuous effort}
  \item{markocc}{integer vector distinguishing marking occasions (1) from 
      sighting occasions (0)}
      \item{trapID}{character column containing detector names (see Details)}
  \item{\dots}{ other arguments to pass to \code{read.table} }
}

\details{

Reads a text file in which the first column is a character string (see Note)
identifying a detector and the next two columns are its x- and
y-coordinates, separated by white space. The coordinates optionally may
be followed by a string of codes `0' or `1' indicating whether the
detector was operated on each occasion. Trap-specific covariates may be
added at the end of the line preceded by `/'. This format is compatible
with the Density software (Efford 2012), except that all detectors are
assumed to be of the same type (usage codes greater than 1 are treated
as 1), and more than one covariate may be specified.

If \code{file} is missing then x-y coordinates will be taken instead
from \code{data}, which should include columns `x' and `y'. Row names of 
\code{data} are read as detector identifiers unless \code{trapID} is specified.
This option does not allow for \code{covariates} or \code{usage}, 
but they maybe added later.

\code{detector} specifies the behaviour of the detector following Efford
et al. (2009). `single' refers to a trap that is able to catch at most
one animal at a time; `multi' refers to a trap that may catch more than
one animal at a time. For both `single' and `multi' detectors a trapped
animals can appear at only one detector per occasion. Detectors of type
`proximity', such as camera traps and hair snags for DNA sampling, allow
animals to be recorded at several detectors on one occasion. See
\code{\link{detector}} for further detector types.

For polygon and transect detector types, each line corresponds to a
vertex and starts with a code to identify the polygon or transect (hence
the same code appears on 2 or more lines). For input from a dataframe
the code column should be named `polyID'. Also, usage and covariates
are for the polygon or transect as a whole and not for each
vertex. Usage and covariates are appended to the end of the line, just
as for point detectors (traps etc.). The usage and covariates for each
polygon or transect are taken from its first vertex. Although the
end-of-line strings of other vertices are not used, they cannot be blank
and should use the same spacing as the first vertex.

}

\value{

An object of class \code{traps} comprising a data frame of x- and
y-coordinates, the detector type (`single', `multi', `proximity',
`count', `polygon' etc.), and possibly other attributes.

}

\references{ 

Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially explicit
capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand. \url{https://www.otago.ac.nz/density/}.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

}

\note{

  Detector names, which become row names in the traps object, should not
  contain underscores.
  
  Prior to 4.3.1 the function did not read usage or covariates from xls 
  or data input.

}

\seealso{  \code{\link{traps}}, \code{\link{make.grid}}, \code{\link{detector}} }

\examples{

\dontrun{
## "trap.txt" should have lines like this 
# 1     365     365
# 2     365     395
# 3     365     425
# etc.
## in following, replace file name with a valid local name
filename <- paste0(system.file("extdata", package = "secr"), '/trap.txt')
tr1 <- read.traps (filename, detector = "proximity")
summary(tr1)

## Or if we have a dataframe of coordinates...
mytrapdf <- data.frame(x = c(365,365,365), y = c(365,395,425),
                       row.names = c('A','B','C'))
mytrapdf
#     x   y
# A 365 365
# B 365 395
# C 365 425
## ...then we can convert it to a `traps' object with
tr2 <- read.traps(data = mytrapdf)
summary(tr2)

}

}

