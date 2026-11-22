\name{verify}
\alias{verify}
\alias{verify.default}
\alias{verify.traps}
\alias{verify.capthist}
\alias{verify.mask}

\title{
Check SECR Data
}

\description{
Check that the data and attributes of an object are internally consistent to avoid crashing functions such as \code{secr.fit}
}

\usage{
\method{verify}{default}(object, report, ...)
\method{verify}{traps}(object, report = 2, ...)
\method{verify}{capthist}(object, report = 2, tol = 0.01, ...)
\method{verify}{mask}(object, report = 2, ...)
}

\arguments{
  \item{object}{
an object of class `traps', `capthist' or `mask'
}
  \item{report}{
integer code for level of reporting to the console. 0 = no report, 1 = errors only, 2 = full.
}
  \item{tol}{
numeric tolerance for deviations from transect line (m)
}
  \item{...}{ 
other arguments (not used) 
}

}
\details{

Checks are performed specific to the class of `object'. The default
method is called when no specific method is available (i.e. class not
`traps', `capthist' or `mask'), and does not perform any checks.

\code{verify.capthist}
\enumerate{
  \item No `traps' component
  \item Invalid `traps' component reported by verify.traps
  \item No live detections (not required for unmarked data, or for telemetry type `marking' with
  sightings and no animals)
  \item Missing values not allowed in capthist
  \item Signal detector signal(s) less than threshold or invalid threshold; `signal' attribute
  does not match detections
  \item Live detection(s) after reported dead
  \item Empty detection histories (except concurrent or marking telemetry and all-sighting data)
  \item More than one detection per detector per occasion at binary detector(s) (`single',
  `multi', `proximity', `polygonX', `transectX', `signal' or `signalnoise')
  \item More than one capture in single-catch trap(s)
  \item Animal detected at more than one detector on an occasion (`multi', `polygonX' or
  `transectX' detectors)
  \item Counts less than zero (`count', `polygon' or `transect' detectors)
  \item Number of rows in `traps' object not compatible with reported detections
  \item Number of rows in dataframe of individual covariates differs from capthist
  \item Number of occasions in usage matrix differs from capthist
  \item Detections at unused detectors
  \item Number of coordinates does not match number of detections
  (`polygon', `polygonX', `transect' or `transectX' detectors)
  \item Coordinates of detection(s) outside polygons (`polygon' or `polygonX' detectors)
  \item Coordinates of detection(s) do not lie on any transect
  (`transect' or `transectX' detectors)
  \item Row names (animal identifiers) not unique, or missing
  \item Inconsistencies in resighting data: dimensions of `markocc', `Tu' or `Tm' incompatible
  with the capthist, or counts negative or not whole numbers; sightings at unused detectors;
  new animals on a sighting-only occasion; sightings on a marking-only occasion
  \item Inconsistencies in telemetry data: occasions with no detections; number or names of
  telemetered animals, or number of locations, not matching the detections; for `dependent'
  telemetry, animals with an all-zero detection history; for `marking' telemetry, a detection
  history without telemetry locations or a telemetered animal without a detection history; for
  `independent' telemetry, non-zero detection histories
  \item Invalid non-target data
}

For multi-session data each session is checked separately, and the following checks are also
made across sessions: levels of factor covariate(s) differ between sessions (a warning, not an
error; see \code{\link{shareFactorLevels}}); individual covariates present in some sessions and
NULL in others (stops with an error).

\code{verify.traps}
\enumerate{
  \item Missing detector coordinates not allowed
  \item Number of rows in dataframe of detector covariates differs from expected
  \item Invalid combination of detector types (a mixture of types is allowed only for point
  detectors -- `single', `multi', `proximity', `count', `signal', `signalnoise', `unmarked'
  or `presence' -- with or without `telemetry')
  \item Number of detectors in usage matrix differs from expected
  \item Number of occasions in usage matrix differs from the length of the detector attribute
  \item Occasions with no used detectors
  \item Number of occasions in \code{markocc} differs from usage matrix (resighting data)
  \item Values in \code{markocc} other than -2, -1, 0 or 1 (resighting data)
  \item Polygons overlap, or no search area specified (`polygon' or `polygonX' detectors)
  \item PolyID missing or not factor (`polygon' or `polygonX' detectors)
  \item Polygon concave in east-west direction (`polygon' or `polygonX' detectors)
  \item More than one row in the notional traps object for telemetry (telemetry-only data)
}

For multi-session data each session is checked separately, and the following checks are also
made across sessions: levels of factor trap covariate(s) differ between sessions (a warning);
trap covariates present in some sessions and NULL in others (stops with an error).

\code{verify.mask}
\enumerate{
  \item Valid x and y coordinates
  \item Number of rows in covariates dataframe differs from expected
}

For multi-session data each session is checked separately, and levels of factor mask
covariate(s) are also compared between sessions (a warning).

Earlier errors may mask later errors: fix & re-run.

}

\value{ 

A list with the component \code{errors}, a logical value indicating
whether any errors were found. If \code{object} contains multi-session
data then session-specific results are contained in a further list
component \code{bysession}. For a single-session capthist the list also has components
\code{trapcheck} (the result of \code{verify.traps}) and \code{captcheck} (a list of logical
values, one for each group of checks on the detections), and, if there are detections at unused
detectors, \code{detections.at.unused.detectors}. For a traps object the components include the
logical result of each check. The value is returned invisibly.

Full reporting is the same as `errors only' except that a message is
posted when no errors are found.

}

\seealso{

\code{\link{capthist}}, \code{\link{make.marking}}, \code{\link{secr.fit}},
\code{\link{shareFactorLevels}}

}

\examples{

verify(captdata)

## create null (complete) usage matrix, and mess it up
temptraps <- make.grid()
usage(temptraps) <- matrix(1, nr = nrow(temptraps), nc = 5)
usage(temptraps)[,5] <- 0
verify (temptraps)

## create mask, and mess it up
tempmask <- make.mask(temptraps)
verify(tempmask)
tempmask[1,1] <- NA
verify(tempmask)

}

