\name{RMarkInput}
\alias{RMarkInput}
\alias{unRMarkInput}

\title{ Convert Data to RMark Input Format }

\description{
  
A single-session capthist object is formed by RMarkInput into a dataframe that may be
passed directly to RMark. 

}

\usage{

RMarkInput(object, grouped = FALSE, covariates = TRUE)

unRMarkInput(df, covariates = TRUE)

}

\arguments{

  \item{object}{secr capthist object}

  \item{grouped}{logical for whether to replace each group of identical
    capture histories with a single line}

  \item{covariates}{logical or character vector; see Details}
  
  \item{df}{dataframe with fields `ch' and `freq' }

}
\details{

  To convert a multi-session object first collapse the sessions with
  \code{\link{join}}.

  If \code{covariates} is TRUE then all columns of individual covariates
  in the input are appended as columns in the output. If
  \code{covariates} is a character-valued vector then only the specified
  covariates will be appended.

  If both \code{grouped} and \code{covariates} are specified in
  \code{RMarkInput}, \code{grouped} will be ignored, with a warning.

}

\value{

  For RMarkInput --
  
  Dataframe with fields ch and freq. `ch' is a character string of 0's
  and 1's. If \code{grouped = FALSE} the rownames are retained and the
  value of `freq' is 1 or -1. Negative values of `freq' indicate
  removal.

  The dataframe also includes individual covariates specified with
  \code{covariates}.

  The attribute `intervals' is copied from `object', if present;
  otherwise it is set to a vector of zeros (indicating a
  closed-population sample).

  For unRMarkInput --

  A single-session capthist object with no traps attribute and hence no
  detector type (i.e. non-spatial capture histories). Covariates are
  copied as requested. 
  
  From \pkg{secr} 4.6.6, missing values (.) in input capture histories are 
  converted to NA in the output, with a warning. The resulting capthist is 
  unusable until the NAs are removed.
  
}

\references{

  Laake, J. and Rexstad E. (2008) Appendix C. RMark - an alternative
approach to building linear models in MARK. In: Cooch, E. and White,
G. (eds) Program MARK: A Gentle Introduction. 6th edition. Available at
\url{http://www.phidot.org/software/mark/docs/book/}.

}

\seealso{
\code{\link{join}}

}

\note{

  In versions before 2.4.0, a spurious occasion was added by
  \code{RMarkInput} when grouped = FALSE. Thanks to Jeff Stetz for
  spotting this.

  The default value for \code{grouped} changed to FALSE in secr 2.4.0.

}

\examples{

## ovenCH is a 5-year mist-netting dataset
ovenRD <- RMarkInput (join(ovenCH))
head(ovenRD)

unRMarkInput(ovenRD)

RMarkInput(deermouse.ESG, covariates = FALSE, grouped = TRUE)
RMarkInput(deermouse.ESG, covariates = TRUE)

\dontrun{
## fit robust-design model in RMark (MARK must be  installed)
library(RMark)
MarkPath <- 'c:/MARK'    ## adjust for your installation
ovenRD.data <- process.data(ovenRD, model = "Robust",
    time.interval = attr(ovenRD, "intervals"))
ovenRD.model <- mark(data = ovenRD.data, model = "Robust",
    model.parameters = list(p = list(formula = ~1, share = TRUE),
    GammaDoublePrime = list(formula = ~1),
    GammaPrime = list(formula = ~1),
    f0 = list(formula = ~1)))   
cleanup(ask = FALSE)
}


}

\keyword{ manip }
