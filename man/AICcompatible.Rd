\name{AICcompatible}
\alias{AICcompatible.secr}
\alias{AICcompatible.secrlist}
\alias{AICcompatible}

\title{Model Compatibility}

\description{

Determine whether models can be compared by AIC. Incompatibility may be due to difference in the data or the specifications of the groups, hcov or binomN arguments to \code{\link{secr.fit}}, 

}

\usage{

\method{AICcompatible}{secr}(object, ...)
\method{AICcompatible}{secrlist}(object, ...)

}

\arguments{

  \item{object}{\code{secr} object output from the function
      \code{\link{secr.fit}}, or a list of such objects with class c("secrlist", "list")}
      
  \item{\dots}{other \code{secr} objects}
  
}

\details{

The capthist objects are checked for strict identity with the function \code{\link{identical}}.

All elements in the output must be TRUE for valid AIC comparison or model averaging using AIC or AICc.

}

\value{

Named logical vector with elements `data', `CL', `groups', `hcov' and `binomN'. 

}

\seealso{ 

\code{\link{AIC.secr}}, \code{\link{modelAverage}}

}

\examples{

AICcompatible(secrdemo.0, secrdemo.CL)

\dontrun{

## A common application of AICcompatible() is to determine 
## the compatibility of models fitted with and without the 
## fastproximity option.

ovenCHp1 <- reduce(ovenCHp, by = 'all', outputdetector = 'count')
ob1 <- secr.fit(ovenCHp, buffer = 300, details = list(fastproximity = TRUE))
ob2 <- secr.fit(ovenCHp1, buffer = 300, details = list(fastproximity = FALSE))
ob3 <- secr.fit(ovenCHp1, buffer = 300, details = list(fastproximity = FALSE), binomN = 1)
AICcompatible(ob1,ob2)
AICcompatible(ob1,ob3)

}

}

\keyword{ models }
