\name{Miscellaneous}

\alias{add.cl}
\alias{detectorcode}
\alias{Dfn2}
\alias{gethazard}
\alias{reparameterize}
\alias{stdform}
\alias{prepareSessionData}

\title{Miscellaneous Internal Functions}

\description{

These functions are exported for possible use elsewhere. They include wrappers for 
C++ code that may be useful in other packages. They are not intended to be 
called by most users.

}

\usage{

add.cl(df, alpha, loginterval, lowerbound = 0) 

detectorcode (object, MLonly = TRUE, noccasions = NULL)

Dfn2(designD, beta = NULL, ...)

gethazard (m, binomN, cc, hk, PIA, usge)

prepareSessionData (capthist, mask, maskusage, design, design0, detectfn, 
    groups, fixed, hcov, details, aslist = TRUE, sessnum = 1) 

reparameterize (realparval, detectfn, details, mask, traps, D, s)

stdform (flist) 

}

\arguments{

add.cl()
  
  \item{df}{data.frame of estimates}
  \item{alpha}{numeric alpha level for confidence intervals}
  \item{loginterval}{logical; if TRUE the interval is a Wald interval on the log scale}
  \item{lowerbound}{numeric lower sutoff for interval }

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
  
prepareSessionData()

\item{capthist}{secr capthist object}
\item{mask}{ \code{\link{mask}} object}
\item{maskusage}{matrix n x m, or NULL}
\item{design}{list}
\item{design0}{list}
\item{detectfn}{integer code for detection function}
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
  
}

\details{


Not for general use. Refer to secr code.

\code{reparameterize} applies the alternative parameterisations in 
\url{https://www.otago.ac.nz/density/pdfs/secr-parameterisations.pdf}.

  \code{Dfn2} is supplied automatically as 'details' argument Dfn in 
  \code{\link{secr.fit}} when the switch Dlambda is set to TRUE for the 
  multi-session trend reparameterization of density. \code{Dfn2} 
  uses beta = NULL to return the required number of density coefficients 
  (beta parameters) in the model.
}

\value{

add.cl -- dataframe with additional columns for lower and upper limits

detectorcode -- integer code for detector type

\code{Dfn2} -- Vector of density values on the link scale, suitable for the internal array (mask x groups x sessions).

gethazard --

prepareSessionData -- list of data required by likelihood fn

reparameterize -- matrix of recast parameter values for alternative parameterisation

stdform -- list of formulae in standard form (named, no LHS)

}

\seealso{

\code{\link{predictDlambda}},
\href{https://www.otago.ac.nz/density/pdfs/secr-trend.pdf}{secr-trend.pdf},

}


