\name{LR.test}
\alias{LR.test}
\title{ Likelihood Ratio Test }

\description{
  Compute likelihood ratio test to compare two fitted models, one nested within the other.
}
\usage{
LR.test(model1, model2)
}
\arguments{
  \item{model1}{ fitted model }
  \item{model2}{ fitted model }
}
\details{

  The fitted models must be of a class for which there is a logLik
  method (e.g., `secr' or `lm'). Check with \code{methods("logLik")}.
  
  The models must be nested (no check is performed - this is up to the
  user), but either model1 or model2 may be the more general model.

  The models must also be compatible by the criteria of \code{\link{AICcompatible}}. 

  The test statistic is twice the difference of the maximized
  likelihoods. It is compared to a chi-square distribution with df equal
  to the number of extra parameters in the more complex model. 

}

\value{
  Object of class `htest', a list with components
  \item{statistic }{value the test statistic }
  \item{parameter }{degrees of freedom of the approximate chi-squared distribution of the test statistic }
  \item{p.value }{probability of test statistic assuming chi-square distribution}
  \item{method }{character string indicating the type of test performed }
  \item{data.name }{character string with names of models compared }
}

\seealso{ 

\code{\link{AICcompatible}},
\code{\link{AIC.secr}}, 
\code{\link{score.test}}  

}

\examples{

## two pre-fitted models
AIC (secrdemo.0, secrdemo.b)
LR.test  (secrdemo.0, secrdemo.b)

}
