\name{read.telemetry}
\alias{read.telemetry}

\title{ Import Telemetry Fixes }

\description{

  A shortcut function for constructing a telemetry capthist object from a file of
telemetry fixes. Telemetry data are generally similar in format to polygon data
(see also \code{\link{addTelemetry}}).

}

\usage{
read.telemetry(file = NULL, data = NULL, covnames = NULL, verify = TRUE, ...)
}

\arguments{
  \item{file}{character name of text file }
  \item{data}{data.frame containing coordinate data (alternative to \code{file}) }
  \item{covnames}{character vector of names for individual covariates }
  \item{verify}{logical for whether to check input }
  \item{\dots}{other arguments passed to countfields, read.table etc. }
}

\details{

Input data may be in a text file (argument \code{file}) or a dataframe
(argument \code{data}). Data should be in the XY format for function `read.capthist` 
i.e. the first 5 columns should be Session, ID, Occasion, X, Y. Further columns are 
treated as individual covariates.

No `traps' input is required. A traps object is generated automatically.

}

\value{
  
An secr capthist object including attribute `telemetryxy' with the x-y
coordinates, and a `traps' object with detector type = `telemetry'

}

\seealso{
\code{\link{addTelemetry}}, \code{\link{read.capthist}}
}

\examples{

\donttest{

## example file of telemetry data included in secr
filename <- system.file("extdata/telemetrydemo.txt", package = "secr")

## peek at raw data
head(readLines(filename))

teCH <- read.telemetry(filename)

plot(teCH)
head(teCH)
RPSV(teCH, CC = TRUE)

fit <- secr.fit(teCH, detectfn = 'HHN', trace = FALSE,
    details = list(safeLL = TRUE, uselog = TRUE), buffer = 100)
predict(fit)

}

}

