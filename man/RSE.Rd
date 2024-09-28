\name{RSE}
\alias{RSE}

\title{RSE from Fitted Model}

\description{

Precision of parameter estimates from an SECR model, expressed as relative standard error.

}

\usage{

RSE(fit, parm = NULL, newdata = NULL)

}

\arguments{

  \item{fit}{secr or openCR fitted model}
  \item{parm}{character; names of one or more real parameters (default all)}
  \item{newdata}{dataframe of covariates for \code{\link{predict.secr}}}

}

\details{

The relative standard error (RSE) of parameter \eqn{\theta} is \eqn{RSE(\hat \theta) = \widehat{SE} (\theta) / {\hat \theta}}. 

For a parameter estimated using a log link with single coefficient \eqn{\beta}, the RSE is also \eqn{\mbox{RSE}(\hat \theta) = \sqrt {\exp( \mbox{var}(\beta))-1}}{RSE(\hat \theta) = \sqrt {exp(var (\beta))-1}}.
This formula is used wherever applicable.

}

\value{

Named vector of RSE, or matrix if newdata has more than one row.

}

\references{

Efford, M. G. and Boulanger, J. 2019. Fast evaluation of study designs for spatially explicit capture--recapture. \emph{Methods in Ecology and Evolution} \bold{10}, 1529--1535.

}

\note{

The less explicit abbreviation CV has been used for the same quantity (sometimes expressed as a percentage). CV is used also for the relative standard deviation of a distribution.

}

\seealso{

\code{\link{CV}}

}

\examples{

RSE(secrdemo.0)

}

