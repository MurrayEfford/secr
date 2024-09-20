\name{Miscellaneous}

\alias{add.cl}
\alias{gethazard}
\alias{reparameterize}
\alias{stdform}
\alias{detectorcode}

\title{Miscellaneous Internal Functions}

\description{

These functions are exported for possible use elsewhere. They include wrappers for 
C++ code that may be useful in other packages. They are not intended to be 
called by most users.

}

\usage{

add.cl(df, alpha, loginterval, lowerbound = 0) 

detectorcode (object, MLonly = TRUE, noccasions = NULL)

gethazard (m, binomN, cc, hk, PIA, usge)

reparameterize (realparval, detectfn, details, mask, traps, D, s)

stdform (flist) 

}

\arguments{

  \item{m}{integer number of mask rows}
  \item{binomN}{integer vector indicating detector type (-2 = multi-catch trap)}
  \item{cc}{integer number of parameter combinations}
  \item{hk}{numeric pre-computed hazards}
  \item{PIA}{parameter index array}
  \item{usge}{numeric matrix of detector usage}
  
  \item{realparval}{matrix of parameter values, one row per combination}
  \item{detectfn}{integer code for detection function}
  \item{details}{list of additional settings, including code for new paraeterisation}
  \item{mask}{ \code{\link{mask}} object}
  \item{traps}{\code{\link{traps}} object}
  \item{D}{numeric scalar density}
  \item{s}{integer occasion}

  \item{flist}{ list of formulae }
  
  \item{object}{traps object}
  \item{MLonly}{logical; if TRUE restrict to ML }
  \item{noccasions}{integer number of occasions }
  
  \item{df}{data.frame of estimates}
  \item{alpha}{numeric alpha level for confidence intervals}
  \item{loginterval}{logical; if TRUE the interval is a Wald interval on the log scale}
  \item{lowerbound}{numeric lower sutoff for interval }
}

\details{

\code{reparameterize} applies the alternative parameterisations in 
\url{https://www.otago.ac.nz/density/pdfs/secr-parameterisations.pdf}.

Not for general use. Refer to secr code.

}

\value{

add.cl -- dataframe with additional columns for lower and upper limits

detectorcode -- integer code for detector type

gethazard --

reparameterize -- matrix of recast parameter values for alternative parameterisation

stdform -- list of formulae in standard form (named, no LHS)


}

