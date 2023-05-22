\name{secr-package}
\alias{secr-package}
\alias{secr}

\docType{package}

\title{Spatially Explicit Capture--Recapture Models }

\description{

  Functions to estimate the density and size of a spatially
  distributed animal population sampled with an array of passive
  detectors, such as traps, or by searching polygons or transects.

}

\details{
\tabular{ll}{
Package: \tab secr\cr
Type: \tab Package\cr
Version: \tab 4.6.1\cr
Date: \tab 2023-05-23\cr
License: \tab GNU General Public License Version 2 or later\cr
}

Spatially explicit capture--recapture is a set of methods for studying
marked animals distributed in space. Data comprise the locations of
detectors (traps, searched areas, etc. described in an object of class
`traps'), and the detection histories of individually marked
animals. Individual histories are stored in an object of class
`capthist' that includes the relevant `traps' object.

Models for population density (animals per hectare) and detection are
defined in \pkg{secr} using symbolic formula notation. Density models
may include spatial or temporal trend. Possible predictors for detection
probability include both pre-defined variables (t, b, etc.)
corresponding to `time', `behaviour' and other effects), and
user-defined covariates of several kinds. Habitat is distinguished from
nonhabitat with an object of class `mask'.

Models are fitted in \pkg{secr} by maximizing either the full likelihood
or the likelihood conditional on the number of individuals observed
(\eqn{n}). Conditional likelihood models are limited to homogeneous
Poisson density, but allow continuous individual covariates for
detection. A model fitted with \code{\link{secr.fit}} is an object
of class \code{secr}. Generic methods (plot, print, summary, etc.) are
provided for each object class.

A link at the bottom of each help page takes you to the
help index. Several vignettes complement the help pages:

\tabular{ll}{
General interest \tab \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-overview.pdf}{secr-overview.pdf} \tab general introduction \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-datainput.pdf}{secr-datainput.pdf} \tab data formats and input functions \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-version4.pdf}{secr-version4.pdf} \tab changes in secr 4.0 \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-manual.pdf}{secr-manual.pdf} \tab consolidated help pages \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-tutorial.pdf}{secr-tutorial.pdf} \tab introductory tutorial \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-habitatmasks.pdf}{secr-habitatmasks.pdf} \tab buffers and habitat masks \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-spatialdata.pdf}{secr-spatialdata.pdf} \tab using spatial data \cr   
  \href{https://www.otago.ac.nz/density/pdfs/secr-models.pdf}{secr-models.pdf} \tab linear models in \pkg{secr} \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-troubleshooting.pdf}{secr-troubleshooting.pdf}  \tab problems with secr.fit, including speed issues \cr
\cr
\cr
More specialised topics \tab \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-densitysurfaces.pdf}{secr-densitysurfaces.pdf} \tab modelling density surfaces \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-finitemixtures.pdf}{secr-finitemixtures.pdf} \tab mixture models for individual heterogeneity \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-markresight.pdf}{secr-markresight.pdf} \tab mark-resight data and models \cr
    \href{https://www.otago.ac.nz/density/pdfs/secr-multisession.pdf}{secr-multisession.pdf} \tab multi-session capthist objects and models \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-noneuclidean.pdf}{secr-noneuclidean.pdf} \tab non-Euclidean distances \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-parameterisations.pdf}{secr-parameterisations.pdf} \tab alternative parameterisations sigmak, a0 \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-polygondetectors.pdf}{secr-polygondetectors.pdf} \tab using polygon and transect detector types \cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-sound.pdf}{secr-sound.pdf} \tab analysing data from microphone arrays\cr
  \href{https://www.otago.ac.nz/density/pdfs/secr-varyingeffort.pdf}{secr-varyingeffort.pdf} \tab variable effort in SECR models\cr
}

The datasets \code{\link{captdata}} and \code{\link{ovenbird}} include examples of fitted
models. For models fitted to other datasets see \href{https://www.otago.ac.nz/density/pdfs/secr-version4.pdf}{secr-version4.pdf} Appendix 2.

Add-on packages extend the capability of \pkg{secr} and are
documented separately. 
\pkg{secrlinear} enables the estimation of linear
density (e.g., animals per km) for populations in linear habitats such
as stream networks
(\href{https://cran.r-project.org/package=secrlinear/vignettes/secrlinear-vignette.pdf}{secrlinear-vignette.pdf}). 
\pkg{secrdesign}
enables the assessment of alternative study designs by Monte Carlo
simulation; scenarios may differ in detector (trap) layout, sampling
intensity, and other characteristics
(\href{https://www.otago.ac.nz/density/pdfs/secrdesign-vignette.pdf}{secrdesign-vignette.pdf}).
\pkg{ipsecr} fits some awkward models (e.g., for single-catch traps) by simulation and inverse prediction 
(\href{https://CRAN.R-project.org/package=ipsecr/vignettes/ipsecr-vignette.pdf}{ipsecr-vignette.pdf}).
\pkg{openCR} fits open population models, both non-spatial and spatial
(\href{https://CRAN.R-project.org/package=openCR/vignettes/openCR-vignette.pdf}{openCR-vignette.pdf}).
 
The analyses in \pkg{secr} extend those available in the software
Density (see
\href{https://www.otago.ac.nz/density/}{www.otago.ac.nz/density/}
for the most recent version of Density).  Help is available on the
`DENSITY | secr' forum at
\href{http://www.phidot.org/forum/}{www.phidot.org} and
the Google group \href{https://groups.google.com/forum/#!forum/secrgroup}{secrgroup}. Feedback
on the software is also welcome, including suggestions for additional
documentation or new features consistent with the overall design.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Borchers, D. L. and Fewster, R. M. (2016)  Spatial capture--recapture models.  
\emph{Statistical Science} \bold{31}, 219--232.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area searches. \emph{Ecology}
\bold{92}, 2202--2207.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

Efford, M. G., Borchers D. L. and Mowat, G. (2013) Varying effort in
  capture--recapture studies. \emph{Methods in Ecology and Evolution}
  \bold{4}, 629--636.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
\emph{Animal Biodiversity and Conservation} \bold{27},
 217--228.

Efford, M. G. and Fewster, R. M. (2013) Estimating population
  size by spatially explicit capture--recapture. \emph{Oikos}
  \bold{122}, 918--928.
  
Efford, M. G. and Hunter, C. M. (2017) Spatial capture--mark--resight
estimation of animal population density. \emph{Biometrics} \bold{74}, 411--420.

Efford, M. G. and Mowat, G. (2014) Compensatory heterogeneity in
capture--recapture data.\emph{Ecology} \bold{95}, 1341--1348.

Royle, J. A., Chandler, R. B., Sollmann, R. and Gardner, B. (2014)
\emph{Spatial capture--recapture.} Academic Press.

Royle, J. A. and Gardner, B. (2011) Hierarchical spatial
capture--recapture models for estimating density from trapping
arrays. In: A.F. O'Connell, J.D. Nichols and K.U. Karanth (eds)
\emph{Camera Traps in Animal Ecology: Methods and Analyses}. Springer,
Tokyo.  Pp. 163--190.
 
}

\keyword{ package }

\seealso{
\code{\link{read.capthist}},
\code{\link[secr]{secr.fit}},
\code{\link[secr]{traps}},
\code{\link[secr]{capthist}},
\code{\link[secr]{mask}}
}
\examples{

\dontrun{

## generate some data & plot
detectors  <- make.grid (nx = 10, ny = 10, spacing = 20,
    detector = "multi")
plot(detectors, label = TRUE, border = 0, gridspace = 20)
detections <- sim.capthist (detectors, noccasions = 5,
    popn = list(D = 5, buffer = 100),
    detectpar = list(g0 = 0.2, sigma = 25))
session(detections) <- "Simulated data"
plot(detections, border = 20, tracks = TRUE, varycol = TRUE)

## generate habitat mask
mask <- make.mask (detectors, buffer = 100, nx = 48)

## fit model and display results
secr.model <- secr.fit (detections, model = g0~b, mask = mask)
secr.model

}

}
