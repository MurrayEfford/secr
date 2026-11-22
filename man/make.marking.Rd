\name{make.marking}
\alias{make.marking}

\title{Capture History for Telemetry-calibrated Cue Counts}

\description{

Build in one step a \code{capthist} object for telemetry type 'marking' (telemetry-calibrated cue
counts; see \code{\link{addTelemetry}}) from the detections and telemetry of collared animals and
the sightings of animals that are not individually identified. The function may be used for
single-session or multi-session data, and allows a session in which no animal was collared.

}

\usage{

make.marking(traps, detections = NULL, telemetry = NULL, marked = NULL,
    covariates = NULL, Tu = NULL, Tm = NULL, noccasions = NULL, verify = TRUE, ...)

}

\arguments{

  \item{traps}{\code{traps} object for the cameras, without a telemetry detector, with
    \code{\link{usage}} (K x S) if effort varies; or, for a multi-session capthist, a named list of
    \code{traps} objects (class c("traps", "list")) one per session. The detector type is taken
    from \code{traps}: "count", or "proximity" (see Details).}

  \item{detections}{data frame of the detections of identified (collared) animals at cameras in the
    format of the \code{captures} argument of \code{\link{make.capthist}} with
    \code{fmt = "trapID"}: columns Session, ID, Occasion, Detector (a row name of \code{traps}),
    one row per detection (repeat rows for counts). Session is the first column, also for a
    single session. May be NULL or have no rows for a session. Column names are ignored; further
    columns are dropped. Only this `trapID' format is accepted (the `XY' format is for polygon
    detectors).}

  \item{telemetry}{telemetry fixes: a data frame with columns ID, x, y (and Session for several
    sessions; a column Occ is ignored); the output of \code{\link{read.telemetry}}; or a list
    of such data frames, one per session. May be NULL, or have no fixes for a session, in which
    case that session has no collared animals.}

  \item{marked}{matrix of the proportion of each sighting occasion that each collared animal was
    detectable (collar active), with one row per collared animal (row names are IDs) and S columns;
    a list of matrices for several sessions (NULL allowed). If NULL all collars are assumed active
    throughout (see \code{\link{marked}}).}

  \item{covariates}{data frame of individual covariates of the collared animals, with row names
    equal to IDs or with a column ID (and Session for several sessions); a list of data frames for
    several sessions. A session without collars receives a data frame with no rows and the same
    columns.}

  \item{Tu, Tm}{sightings of unmarked animals (\code{Tu}) and of marked animals that were not
    identified (\code{Tm}) as counts for each detector and occasion (K x S); see Details.}

  \item{noccasions}{integer number of sighting occasions S; needed only if \code{traps} has no
    usage attribute and \code{Tu} is not a matrix.}

  \item{verify}{logical; if TRUE the result is checked with \code{\link{verify}}.}

  \item{\dots}{other arguments passed to \code{\link{read.table}} when reading files of sightings.}

}

\details{

The usual route to a capthist for telemetry type 'marking' is a sequence of functions:
\code{\link{make.capthist}} for the detections of collared animals, \code{\link{read.telemetry}},
\code{\link{addTelemetry}} with \code{type = "marking"}, \code{\link{marked<-}},
\code{\link{addSightings}} and \code{covariates<-}. \code{make.marking} combines these
steps and produces an object of the same form: the cameras are detectors 1..K, followed by a
notional telemetry detector (a copy of the first camera) that is used on one additional occasion
that pools all fixes; \code{markocc} is 0 for the S sighting occasions and 1 for the telemetry
occasion; the usage matrix is (K + 1) x (S + 1); the attribute \code{telemetryxy} holds the fixes
of each collared animal; and there is one row for every collared animal (every animal with fixes),
all-zero if it was never detected at a camera.

Unlike the sequence, which returns the detections unchanged when there is no telemetry,
\code{make.marking} builds a complete object for a session in which no animal was collared
(no fixes). This is needed for a multi-session analysis in which detection parameters are shared
among sessions and density is estimated for each, including sessions without collars
(\code{D ~ session}). At least one session must have a collared animal.

\emph{Sessions.} For several sessions \code{traps} must be a named list; the names are the session
names, and the Session values in \code{detections} (and in \code{telemetry} and \code{covariates}
if they have a Session column) must match them. For one session the session name is the single
value of Session in \code{detections} (otherwise in \code{telemetry}, otherwise "1"). A data frame
for several sessions may be replaced by a list with one component per session.

\emph{Sightings.} \code{Tu} and \code{Tm} may be given as

\itemize{
  \item a data frame, or the name of a text file (read with \code{\link{read.table}} and the
  \dots arguments), whose first column is a session label and whose other columns are the
  occasions, with K rows (one per detector, in the order of \code{traps}) for each session stacked
  one above another (the format of \code{\link{addSightings}}). Sessions are taken in order of
  first appearance, not alphabetical order. If the labels equal the session names they are matched
  by name; otherwise they are matched by position, with a message.
  \item a matrix (single session only)
  \item a list with one component per session, each NULL, a matrix, or the name of a file whose
  first column is dropped
}

and the dimensions are checked against \code{traps} and its usage. Uncertain sightings (attribute
\code{Tn}) are not accepted because no sighting model uses them.

\emph{Detector types.} Detectors are "count" (expected count is a hazard; the model fitted by
\code{\link{secr.fit}}) or "proximity" (binary detection). A binary model for telemetry-calibrated
cue counts is not yet available in \pkg{secr}, and \code{secr.fit} stops if given proximity
detectors with telemetry type 'marking'; a message is given when proximity detectors are
specified. With proximity detectors repeated detections of an animal at a detector on one
occasion are set to 1, with a warning. With "count" detectors a warning is given if every count in
\code{Tu} is 0 or 1, as that may indicate presence/absence data for which the count model is not
appropriate.

Animals are in sorted order of ID, as from \code{read.telemetry}, and row names are the IDs
(\code{addTelemetry} pads IDs to a common width). Individual covariates are matched to animals by
ID, not by position, and have the IDs as row names.

}

\value{

A \code{capthist} object, or for several sessions a multi-session \code{capthist}, with
telemetry type 'marking' and, if specified, attributes \code{Tu}, \code{Tm} (see
\code{\link{addSightings}}), \code{marked} and individual covariates. Factor covariates have the
same levels in all sessions (\code{\link{shareFactorLevels}}).

}

\seealso{

\code{\link{addTelemetry}}, \code{\link{addSightings}}, \code{\link{make.capthist}},
\code{\link{read.telemetry}}, \code{\link{markocc}}, \code{\link{marked}},
\code{\link{secr.fit}}, \code{\link{verify}}

}

\examples{

## cameras: 5 x 5 grid, 4 sighting occasions
tr <- make.grid(nx = 5, ny = 5, spacing = 100, detector = "count")
usage(tr) <- matrix(1, nrow(tr), 4)

## fixes for three collared animals (A, B, C)
set.seed(123)
tel <- data.frame(ID = rep(c("A","B","C"), each = 10), x = rnorm(30, 200, 40),
    y = rnorm(30, 200, 40))

## detections of two of them at cameras (Session, ID, Occasion, Detector)
det <- data.frame(Session = 1, ID = c("A","A","C"), Occasion = c(1,3,2),
    Detector = rownames(tr)[c(1,7,12)])

## sightings of unmarked animals (detectors x occasions)
Tu <- matrix(rpois(25*4, 0.5), 25, 4)

ch <- make.marking(tr, det, tel, Tu = Tu)
summary(ch)

## two sessions; no animal was collared in the second
trl <- list(year1 = tr, year2 = tr)
class(trl) <- c("traps", "list")
det$Session <- "year1"
tel$Session <- "year1"
Tu2 <- list(Tu, matrix(rpois(25*4, 0.5), 25, 4))
ch2 <- make.marking(trl, det, tel, Tu = Tu2)
sapply(ch2, nrow)   # number of collared animals in each session

}
