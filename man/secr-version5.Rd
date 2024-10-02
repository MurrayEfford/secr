\name{secr-version5}
\alias{secr-version5}
\alias{version5}

\title{ Changes in \pkg{secr} 5.0 }

\description{

This document explains changes in \pkg{secr} 5.0. Version 5.0 is compatible in most respects with earlier versions, but a few names and one default have been changed without warning. 
See the \href{https://CRAN.R-project.org/package=secr/NEWS}{NEWS} file for a 
complete list of the changes over time.

}

\section{Potential Gotchas}{

\subsection{New generics}{

Several new generic functions are defined, with methods specifically for `secr' 
fitted models (\code{\link{esa}}, \code{\link{fxi}}, \code{\link{fxTotal}}).

}

\subsection{New names}{

Some functions with "." in their name have been renamed to avoid confusion 
with methods for generics. 

Where possible, the old names have been deprecated (tagged with a warning), and will
continue to work for a while.

\tabular{ll}{
Old \tab New \cr
 \code{buffer.contour} \tab \code{\link{bufferContour}} \cr
 \code{esa.plot} \tab \code{\link{esaPlot}} \cr
 \code{fxi.contour} \tab \code{\link{fxiContour}} \cr
 \code{fxi.mode} \tab \code{\link{fxiMode}} \cr
 \code{fxi.secr} \tab \code{\link{fxi}} (generic) \cr
 \code{fx.total} \tab \code{\link{fxTotal}} (generic) \cr
 \code{pdot.contour} \tab \code{\link{pdotContour}} \cr
}

\code{fxi.secr} has been replaced by the generic \code{\link{fxi}}. Thus instead of \cr
\code{fxi.secr(secrdemo.0, i = 1, X = c(365,605))} use \cr
\code{fxi(secrdemo.0, i = 1, X = c(365,605))}.
}

\subsection{New default}{

AIC and related functions now default to criterion = "AIC" instead of 
criterion = "AICc".

Some of us have been uneasy for a long time about blanket use of the AICc 
small-sample adjustment to AIC (Hurvich and Tsai 1989). Royle et al. (2014) 
expressed doubts because the sample size itself is poorly defined. AICc is 
widely used, but AIC may be better for model averaging even when samples 
are small (Turek and Fletcher 2012; Fletcher 2019, p. 60). 
}

}

\section{New features}{

\subsection{New data \code{blackbearCH}}{

\pkg{secr} 5.0 includes a new black bear DNA hair snag dataset from the Great Smoky 
Mountains, Tennessee (thanks to J. Laufenberg, F. van Manen and J. Clark).

}

\subsection{Goodness-of-fit \code{MCgof}}{

The method of Choo et al. (2024) for emulating the Bayesian p-value 
goodness-of-fit test (Gelman 1996, Royle et al. 2014) has been implemented 
as the generic \code{\link{MCgof}} with a method for `secr' fitted models. 
I thank Yan Ru Choo for his assistance.

This is a new approach and should be used with caution. Bugs may yet be found, 
and the power of the tests is limited. 

}
}

\section{Extended capability}{

These extensions allow \code{MCgof} to cover a wider range of models:

\itemize{
\item \code{\link{detectpar}} optionally returns values for each detector
\item \code{\link{pdot}} accepts detector- and occasion-specific detection parameters
}

}

\section{Changes behind the scenes}{

The code for area-search and transect-search models (detector types `polygonX',
`polygon', `transectX', `transect') has been streamlined with a view to 
removing it to another package. Simulation for these models (functions
\code{sim.capthist}, \code{sim.detect}) will remain in \pkg{secr}, but uses native R 
functions rather than RcppNumerical of Qiu et al. (2023). 

The undocumented detection function `HPX' has been removed.

}

\seealso{

\link{secr-deprecated}, \link{secr-defunct}

}

\references{

Choo, Y. R., Sutherland, C. and Johnston, A. (2024) \href{https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14386}{A Monte Carlo resampling 
framework for implementing goodness-of-fit tests in spatial capture-recapture models}
\emph{Methods in Ecology and Evolution} DOI: 10.1111/2041-210X.14386.

Efford, M. G. (2024) secr: Spatially explicit capture-recapture models. R package version 5.0.0.  \url{https://CRAN.R-project.org/package=secr}

Fletcher, D. (2019) \emph{Model averaging}. SpringerBriefs in Statistics. Berlin: Springer-Verlag.

Hurvich, C. M. and Tsai, C. L. (1989) Regression and time series model selection in small samples. \emph{Biometrika} \bold{76}, 297--307.

Gelman, A., Meng, X.-L., and Stern, H. (1996) Posterior predictive assessment of 
model fitness via realized discrepancies. \emph{Statistica Sinica} \bold{6}, 733--807.

Qiu, Y., Balan, S., Beall, M., Sauder, M., Okazaki, N. and Hahn, T. (2023) 
RcppNumerical: 'Rcpp' Integration for Numerical Computing Libraries. R package 
version 0.6-0. \url{https://CRAN.R-project.org/package=RcppNumerical}

Royle, J. A., Chandler, R. B., Sollmann, R. and Gardner, B. (2014)
\emph{Spatial capture--recapture}. Academic Press.

Turek, D. and Fletcher, D. (2012) Model-averaged Wald confidence
intervals. \emph{Computational statistics and data analysis} \bold{56},
2809--2815.

}
