\name{collate}
\alias{collate}
\alias{collate.secr}
\alias{collate.secrlist}

\title{ Array of Parameter Estimates }

\description{
  
Estimates from one or more openCR models are formed into an array.

}

\usage{

\method{collate}{secr}(object, ..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, perm = 1:4, fields = 1:4)

\method{collate}{secrlist}(object, ..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, perm = 1:4, fields = 1:4)

}

\arguments{
  \item{object}{secr or secrlist object}
  \item{\dots}{other secr objects } 
  \item{realnames}{ character vector of real parameter names }
  \item{betanames}{ character vector of beta parameter names }
  \item{newdata}{ optional dataframe of values at which to evaluate models }
  \item{alpha}{ alpha level for confidence intervals }
  \item{perm}{ permutation of dimensions in output from \code{collate} }
  \item{fields}{ vector to restrict summary fields in output }
}

\details{

\code{collate} extracts parameter estimates from a set of fitted secr
model objects. \code{fields} may be used to select a subset of summary
fields ("estimate","SE.estimate","lcl","ucl") by name or number. 

}

\value{

A 4-dimensional array of model-specific parameter
estimates. By default, the dimensions correspond respectively to rows in
\code{newdata} (usually sessions), models, statistic fields (estimate, SE.estimate, lcl,
ucl), and parameters ("D", "g0" etc.). For particular comparisons it often helps
to reorder the dimensions with the \code{perm} argument.

}

\seealso{ \code{\link{modelAverage}}, \code{\link{secr.fit}}}

\examples{

collate (secrdemo.0, secrdemo.b, perm = c(4,2,3,1))[,,1,]

}

