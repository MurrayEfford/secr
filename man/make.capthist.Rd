\name{make.capthist}
\alias{make.capthist}

\title{ Construct capthist Object }

\description{
  Form a \code{capthist} object from a data frame of capture records and a \code{traps} object.
}

\usage{
make.capthist(captures, traps, fmt = c("trapID", "XY"), noccasions = NULL,
    covnames = NULL, bysession = TRUE, sortrows = TRUE,
    cutval = NULL, tol = 0.01, snapXY = FALSE, noncapt = "NONE", signalcovariates)
}

\arguments{
  \item{captures}{ dataframe of capture records in one of two possible formats (see Details) }
  \item{traps}{ object of class \code{traps} describing an array of passive detectors }
  \item{fmt}{ character string for capture format. }
  \item{noccasions}{ number of occasions on which detectors were operated }
  \item{covnames}{ character vector of names for individual covariate fields }
  \item{bysession}{ logical, if true then ID are made unique by session }
  \item{sortrows}{ logical, if true then rows are sorted in ascending order of animalID }
  \item{cutval}{ numeric, threshold of signal strength for `signal' detector type }
  \item{tol}{ numeric, snap tolerance in metres}
  \item{snapXY}{ logical; if TRUE then fmt = 'XY' uses nearest trap
  within tol for non-polygon detectors }
  \item{noncapt}{ character value; animal ID used for `no captures' }
  \item{signalcovariates}{character vector of field names from
    `captures'}
}

\details{

\code{make.capthist} is the most flexible way to prepare data for
\code{secr.fit}. See \code{\link{read.capthist}} for a more streamlined
way to read data from text files for common detector types. Each row of
the input data frame \code{captures} represents a detection on one
occasion. The capture data frame may be formed from a text file with
\code{read.table}.

Input formats are based on the Density software (Efford 2012; see also
\href{../doc/secr-datainput.pdf}{secr-datainput.pdf}). If \code{fmt =
"XY"} the required fields are (session, ID, occasion, x, y) in that
order. If \code{fmt = "trapID"} the required fields are (session, ID,
occasion, trap), where \code{trap} is the numeric index of the relevant
detector in \code{traps}. \code{session} and \code{ID} may be
character-, vector- or factor-valued; other required fields are
numeric. Fields are matched by position (column number), \emph{not} by
name. Columns after the required fields are interpreted as individual
covariates that may be continuous (e.g., size) or categorical (e.g.,
age, sex).

If \code{captures} has data from multiple sessions then \code{traps} may
be either a list of \code{traps} objects, one per session, or a single
\code{traps} object that is assumed to apply throughout. Similarly,
\code{noccasions} may be a vector specifying the number of occasions in
each session.

Covariates are assumed constant for each individual; the first
non-missing value is used. The length of \code{covnames} should equal the
number of covariate fields in \code{captures}.

\code{bysession} takes effect when the same individual is detected in
two or more sessions: TRUE results in one capture history per session,
FALSE has the effect of generating a single capture history (this is not
appropriate for the models currently provided in \pkg{secr}).

Deaths are coded as negative values in the occasion field of
\code{captures}. Occasions should be numbered 1, 2, ..., noccasions. By
default, the number of occasions is the maximum value of `occasion' in
\code{captures}.

Signal strengths may be provided in the fifth (fmt = trapID) or sixth
(fmt = XY) columns. Detections with signal strength missing (NA) or
below `cutval' are discarded.

A session may result in no detections. In this case a null line is
included in \code{captures} using the animal ID field given by
\code{noncapt}, the maximum occasion number, and any trapID (e.g. "sess1
NONE 5 1" for a 5-occasion session) (or equivalently "sess1 NONE 5 10
10" for fmt = XY).

Nonspatial data (Session, AnimalID, Occasion and possibly individual 
covariates) may be entered by omitting the `traps' argument or setting 
it to NULL.

}

\value{

An object of class \code{\link{capthist}} (a matrix or array of
detection data with attributes for detector positions etc.). For
`single' and `multi' detectors this is a matrix with one row per animal
and one column per occasion (dim(capthist)=c(nc,noccasions)); each
element is either zero (no detection) or a detector number (the row
number in \code{traps} \emph{not} the row name). For `proximity'
detectors \code{capthist} is an array of values \{-1, 0, 1\} and
dim(capthist)=c(nc,noccasions,ntraps). The number of animals \code{nc}
is determined from the input, as is \code{noccasions} if it is not specified.
\code{traps}, \code{covariates} and other data are retained as
attributes of \code{capthist}.

Deaths during the experiment are represented as negative values in \code{capthist}.

For `signal' and `signalnoise' detectors, the columns of \code{captures}
identified in \code{signalcovariates} are saved along with signal
strength measurements in the attribute `signalframe'.

If the input has data from multiple sessions then the output is an
object of class c("capthist", "list") comprising a list of single-session
\code{capthist} objects.

}

\note{

\code{make.capthist} requires that the data for \code{captures} and
\code{traps} already exist as \R objects. To read data from external
(text) files, first use \code{read.table} and \code{read.traps}, or try
\code{\link{read.capthist}} for a one-step solution.

Prior to \pkg{secr} 4.4.0, occasional valid records for "multi" and "single" detectors 
were rejected as duplicates.

From \pkg{secr} 4.5.0, `snapXY' works for transects as well as point detectors.

}

\references{

Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially explicit
capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand. \url{https://www.otago.ac.nz/density/}.

}

\seealso{ \code{\link{capthist}}, \code{\link{traps}},
\code{\link{read.capthist}}, \code{\link{secr.fit}},
\code{\link{sim.capthist}} }

\examples{

## peek at demonstration data
head(captXY)
head(trapXY)

demotraps <- read.traps(data = trapXY)
demoCHxy  <- make.capthist (captXY, demotraps, fmt = "XY")

demoCHxy            ## print method for capthist
plot(demoCHxy)      ## plot method for capthist
summary(demoCHxy)   ## summary method for capthist


## To enter `count' data without manually repeating rows
## need a frequency vector f, length(f) == nrow(captXY)
n <- nrow(captXY)
f <- sample (1:5, size = n, prob = rep(0.2,5), replace = TRUE)
## repeat rows as required...
captXY <- captXY[rep(1:n, f),]
counttraps <- read.traps(data = trapXY, detector = "count")
countCH  <- make.capthist (captXY, counttraps, fmt = "XY")

}

\keyword{ manip }
