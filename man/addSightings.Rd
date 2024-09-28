\name{addSightings}
\alias{addSightings}

\title{Mark-resight Data}

\description{
  
  Add sighting data on unmarked individuals and/or unidentified marked
  individuals to an existing capthist object. 
  
}

\usage{

addSightings(capthist, unmarked = NULL, nonID = NULL, uncertain = NULL, verify = TRUE,
    ...)

}

\arguments{
  \item{capthist}{secr capthist object}
  
  \item{unmarked}{matrix or list of matrices of sightings of unmarked animals, Tu, or file name (see Details)}
  
  \item{nonID}{matrix or list of matrices of unidentified sightings of marked animals, Tm, or file name (see Details)}

  \item{uncertain}{matrix or list of matrices of uncertain sightings, Tn, or file name (see Details)}
  
  \item{verify}{ logical; if TRUE then the resulting capthist object is
  checked with \code{\link{verify}}}

  \item{\dots}{other arguments passed to \code{\link{read.table}} }

}

\details{

The capthist object for mark-resight analysis comprises distinct marking and sighting occasions, defined in the markocc attribute of \code{traps(capthist)}. Add this attribute to \code{traps(capthist)} with \code{\link{markocc}} before using 'addSightings'. See also \code{\link{read.traps}} and \code{\link{read.capthist}}.

Mark-resight data may be binary (detector type `proximity') or counts (detector types `count', 'polygon' or 'transect'). The detector type is an attribute of \code{traps(capthist)}. Values in \code{unmarked} and \code{nonID} should be whole numbers, and may be greater than 1 even for binary proximity detectors because multiple animals may be detected simultaneously at one place.

Arguments \code{unmarked}, \code{nonID}, \code{uncertain} provide data for attributes
`Tu', `Tm', `Tn' respectively. They may take several forms

\itemize{
  \item a single integer, the sum of all counts*
  \item a matrix of the count on each occasion at each detector (dimensions K x S, where K is the number of detectors and S is the total number of occasions). Columns corresponding to marking occasions should be all-zero.
  \item for multi-session data, a list with components as above
  \item a character value with the name of a text file containing the data; the file will be read with \code{\link{read.table}}. The \dots argument allows some control over how the file is read. The data format comprises at least S+1 columns. The first is a session identifier used to split the file when the data span multiple sessions; it should be constant for a single-session capthist. The remaining S columns contain the counts for occasions 1:S, one row per detector. Further columns may be present; they are ignored at present.    
}

* although this is convenient, the full matrix of counts provides more flexibility (e.g., when you wish to subset by occasion), and enables modelling of variation across detectors and occasions.

}
\value{

  A capthist object with the same structure as the input, but with new sighting-related attributes Tu (sightings of unmarked animals) and/or Tm (unidentified sightings of marked animals). Input values, including NULL, overwrite existing values.

}

\seealso{
    
  \code{\link{markocc}}, 
  \code{\link{read.capthist}}, 
  \code{\link{read.traps}}, 
  \code{\link{sim.resight}},
  \code{\link{Tm}},
  \code{\link{Tu}},
  \code{\link{Tn}},
  \href{https://www.otago.ac.nz/density/pdfs/secr-markresight.pdf}{secr-markresight.pdf} 
  
}

\section{Warning}{** Mark-resight data formats and models are experimental and subject to change **}

\examples{

## construct capthist object MRCH from text files provided in 
## 'extdata' folder, assigning attribute 'markocc' and add unmarked
## and marked sightings from respective textfiles

datadir <- system.file("extdata", package = "secr")
captfile <- paste0(datadir, '/MRCHcapt.txt')
trapfile <- paste0(datadir, '/MRCHtrap.txt')
Tufile <- paste0(datadir, '/Tu.txt')
Tmfile <- paste0(datadir, '/Tm.txt')

MRCH <- read.capthist(captfile, trapfile, detector = c("multi", 
    rep("proximity",4)), markocc = c(1,0,0,0,0))
MRCH1 <- addSightings(MRCH, Tufile, Tmfile)

## alternatively (ignoring marked, not ID sightings)

MRCH <- read.capthist(captfile, trapfile, detector = c("multi", 
    rep("proximity",4)), markocc = c(1,0,0,0,0))
Tu <- read.table(Tufile)[,-1]  # drop session column
MRCH2 <- addSightings(MRCH, unmarked = Tu)
summary(MRCH2)

}

