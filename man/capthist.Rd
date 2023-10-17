\name{capthist} \alias{capthist}

\title{ Spatial Capture History Object }

\description{
A \code{capthist} object encapsulates all data needed by
\code{secr.fit}, except for the optional habitat mask.
}

\details{
An object of class \code{capthist} holds spatial capture histories,
detector (trap) locations, individual covariates and other data needed
for a spatially explicit capture-recapture analysis with
\code{secr.fit}. 

A \code{capthist} is primarily an array of values with dim(capthist) = c(nc,
noccasions, ntraps) where nc is the number of detected individuals. 
Values maybe binary (\{--1, 0, 1\}) or integer depending on the detector type.

Deaths during the experiment are represented as negative values. 

Ancillary data are retained as attributes of a \code{capthist} object as follows:
\itemize{

\item{traps --- object of class \code{traps} (required)}
\item{session --- session identifier (required)}
\item{covariates --- dataframe of individual covariates (optional)}
\item{cutval --- threshold of signal strength for detection (`signal' only)}
\item{signalframe --- signal strength values etc., one row per detection (`signal' only)}
\item{detectedXY --- dataframe of coordinates for location within
  polygon (`polygon'-like detectors only)}
\item{xylist --- coordinates of telemetered animals}
\item{Tu --- detectors x occasions matrix of sightings of unmarked animals}
\item{Tm --- detectors x occasions matrix of sightings of marked but unidentified animals }
\item{Tn --- detectors x occasions matrix of sightings with unknown mark status }

}

\code{read.capthist} is adequate for most data input. Alternatively, the parts of a 
capthist object can be assembled with the function \code{make.capthist}.  
Use \code{sim.capthist} for Monte Carlo simulation
(simple models only). Methods are provided to display and manipulate
\code{capthist} objects (print, summary, plot, rbind, subset, reduce)
and to extract and replace attributes (covariates, traps, xy).

A multi-session \code{capthist} object is a list in which each component
is a \code{capthist} for a single session.  The list maybe derived
directly from multi-session input in Density format, or by combining
existing \code{capthist} objects with \code{\link{MS.capthist}}.

}

\references{ 

Borchers, D. L. and Efford, M. G. (2008) Spatially
explicit maximum likelihood methods for capture--recapture studies.
\emph{Biometrics} \bold{64}, 377--385.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.
}

\note{

Early versions of \pkg{secr} (before 3.0) used an individual x occasion matrix 
for data from single-catch and multi-catch traps, instead of a 3-D array. 
Entries in the matrix corresponded to trap numbers. The function 
\code{\link{updateCH}} converts the old format.

}

\seealso{ \code{\link{traps}}, \code{\link{secr.fit}},
\code{\link{read.capthist}}, \code{\link{make.capthist}},
\code{\link{sim.capthist}},
\code{\link{subset.capthist}}, \code{\link{rbind.capthist}},
\code{\link{MS.capthist}},
\code{\link{reduce.capthist}}, \code{\link{mask}} }

\keyword{ classes }
