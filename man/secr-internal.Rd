\name{Internal}

\alias{add.cl}
\alias{boundarytoSF}
\alias{detectorcode}
\alias{Dfn2}
\alias{gethazard}
\alias{parnames}
\alias{prepareSessionData}
\alias{reparameterize}
\alias{stdform}
\alias{valid.detectfn}

\title{Internal Functions}

\description{ 

Various functions called internally by \pkg{secr} and also exported.

}

\usage{
add.cl(df, alpha, loginterval, lowerbound = 0) 

boundarytoSF (poly)

detectorcode (object, MLonly = TRUE, noccasions = NULL)

Dfn2(designD, beta = NULL, ...)

gethazard (m, binomN, cc, hk, PIA, usge)

parnames (detectfn)

prepareSessionData (capthist, mask, maskusage, design, design0, detectfn, 
    groups, fixed, hcov, details, aslist = TRUE, sessnum = 1) 

reparameterize (realparval, detectfn, details, mask, traps, D, s)

stdform (flist) 

valid.detectfn (detectfn, valid = c(0:3,5:19))

}

\arguments{

add.cl()
  \item{df}{data.frame of estimates}
  \item{alpha}{numeric alpha level for confidence intervals}
  \item{loginterval}{logical; if TRUE the interval is a Wald interval on the log scale}
  \item{lowerbound}{numeric lower sutoff for interval }

boundarytoSF()
  \item{poly}{data defining spatial polygons}

detectorcode()
  \item{object}{traps object}
  \item{MLonly}{logical; if TRUE restrict to ML }
  \item{noccasions}{integer number of occasions }

Dfn2()
  \item{designD}{dataframe of density design data (output from \code{\link{D.designdata}})}
  \item{beta}{numeric vector of beta values (see Details for NULL)}
  \item{\dots}{other arguments (not used)}

gethazard()
  \item{m}{integer number of mask rows}
  \item{binomN}{integer vector indicating detector type (-2 = multi-catch trap)}
  \item{cc}{integer number of parameter combinations}
  \item{hk}{numeric pre-computed hazards}
  \item{PIA}{parameter index array}
  \item{usge}{numeric matrix of detector usage}

parnames()
  \item{detectfn}{integer code for detection function}

prepareSessionData()
  \item{capthist}{secr capthist object}
  \item{mask}{ \code{\link{mask}} object}
  \item{maskusage}{matrix n x m, or NULL}
  \item{design}{list}
  \item{design0}{list}
  \item{groups}{character vector}
  \item{fixed}{list}
  \item{hcov}{character}
  \item{details}{list of additional settings, including code for new paraeterisation}
  \item{aslist}{logical}
  \item{sessnum}{integer}

reparameterize()
  \item{realparval}{matrix of parameter values, one row per combination}
  \item{traps}{\code{\link{traps}} object}
  \item{D}{numeric scalar density}
  \item{s}{integer occasion}

stdform()
  \item{flist}{ list of formulae }
  
valid.detectfn()
  \item{valid}{integer vector of allowed detectfn}

}

\details{

Not for general use. Refer to secr code.

\code{boundarytoSF} coerces various formats to sfc.

\code{reparameterize} applies the alternative parameterisations in 
\url{https://www.otago.ac.nz/density/pdfs/secr-parameterisations.pdf}.

\code{Dfn2} is supplied automatically as 'details' argument Dfn in 
  \code{\link{secr.fit}} when the switch Dlambda is set to TRUE for the 
  multi-session trend reparameterization of density. \code{Dfn2} 
  uses beta = NULL to return the required number of density coefficients 
  (beta parameters) in the model.

}

\value{

\code{add.cl} -- dataframe with additional columns for lower and upper limits

\code{boundarytoSF} -- Spatial object of \pkg{sf} class sfc, containing 
a geometry set of type POLYGON or MULTIPOLYGON. NULL input results in NULL output.

\code{detectorcode} -- integer code for detector type

\code{Dfn2} -- Vector of density values on the link scale, suitable for the internal array (mask x groups x sessions).

\code{gethazard} --

\code{parnames} -- character vector of parameter names

\code{prepareSessionData} -- list of data required by likelihood fn

\code{reparameterize} -- matrix of recast parameter values for alternative parameterisation

\code{stdform} -- list of formulae in standard form (named, no LHS)

\code{valid.detectfn} -- integer code

}

\seealso{

\code{\link{predictDlambda}},
\href{https://www.otago.ac.nz/density/pdfs/secr-trend.pdf}{secr-trend.pdf},

\code{\link{pointsInPolygon}}, 
\href{https://www.otago.ac.nz/density/pdfs/secr-spatialdata.pdf}{secr-spatialdata.pdf}

}

\references{

Hijmans, R. J. (2022) terra: Spatial Data Analysis. R package version 1.5-14. https://rspatial.org/terra/

Pebesma, E. (2018) Simple features for R: standardized support for spatial vector data. \emph{The R Journal} \bold{10(1)}, 439--446.
  https://doi.org/10.32614/RJ-2018-009

Pebesma, E.J. and Bivand, R. S. (2005) Classes and methods for spatial data in R. \emph{R News} \bold{5(2)}, 9--13.
  https://cran.r-project.org/doc/Rnews/Rnews_2005-2.pdf.

}

\examples{

\dontrun{

poly <- cbind(x = c(0,6,6,0,0), y = c(0,0,6,6,0))  
secr::boundarytoSF(poly)

}

}
