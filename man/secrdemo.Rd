\name{secrdemo}
\docType{data}

\alias{secrdemo}
\alias{captXY}
\alias{trapXY}
\alias{captdata}
\alias{secrdemo.0}
\alias{secrdemo.b}
\alias{secrdemo.CL}

\title{ SECR Models Fitted to Demonstration Data }

\description{

  Demonstration data from program Density are provided as text
  files in the `extdata' folder, as raw dataframes (\code{trapXY},
  \code{captXY}), and as a combined \code{capthist} object
  (\code{captdata}) ready for input to \code{secr.fit}.

  The fitted models are objects of class \code{secr} formed by 

      \code{secrdemo.0 <- secr.fit (captdata)}

      \code{secrdemo.b <- secr.fit (captdata, model = list(g0 = ~b))}

      \code{secrdemo.CL <- secr.fit (captdata, CL = TRUE) }

}

\usage{data(secrdemo)}

\details{

The raw data are 235 fictional captures of 76 animals over 5 occasions
in 100 single-catch traps 30 metres apart on a square grid with origin
at (365,365).

Dataframe \code{trapXY} contains the data from the Density input file
`trap.txt', and \code{captXY} contains the data from `capt.txt' (Efford
2012).

The fitted models use a halfnormal detection function and the likelihood
for multi-catch traps (expect estimates of g0 to be biased because of
trap saturation Efford et al.  2009). The first is a null model
(i.e. parameters constant) and the second fits a learned trap response.


\tabular{ll}{
  Object \tab Description \cr
  captXY \tab data.frame of capture data \cr
  trapXY \tab data.frame of trap locations \cr
  captdata \tab capthist object \cr
  secrdemo.0 \tab fitted secr model -- null \cr
  secrdemo.b \tab fitted secr model -- g0 trap response \cr
  secrdemo.CL \tab fitted secr model -- null, conditional likelihood \cr
}
 
}

\references{

Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially explicit
capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand. \url{https://www.otago.ac.nz/density/}.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

}

\seealso{ \code{\link{capthist}}, \code{\link{read.capthist}} }

\examples{

\dontrun{

## navigate to folder with raw data files
olddir <- setwd (system.file("extdata", package="secr"))

## construct capthist object from raw data
captdata  <- read.capthist ("capt.txt", "trap.txt", fmt = "XY", detector = "single")

## generate demonstration fits
secrdemo.0 <- secr.fit (captdata)
secrdemo.CL <- secr.fit (captdata, CL = TRUE)
secrdemo.b <- secr.fit (captdata, model = list(g0 ~ b))

## restore previous setting
setwd(olddir)
}

## display the null model fit, using the print method for secr
secrdemo.0

## compare fit of models
AIC(secrdemo.0, secrdemo.b)

## display estimates for the two models (single session)
collate(secrdemo.0, secrdemo.b)[1,,,]

}

