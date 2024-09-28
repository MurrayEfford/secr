\name{logit}
\alias{logit}
\alias{invlogit}

\title{ Logit Transformation }

\description{
  Transform real values to the logit scale, and the inverse.
}

\usage{
logit(x)
invlogit(y)
}

\arguments{
  \item{x}{ vector of numeric values in (0,1) (possibly a probability) }
  \item{y}{ vector of numeric values }
}

\details{
  The logit transformation is defined as \eqn{\mathrm{logit}(x) = \mathrm{log}( \frac{x}{1-x})}{logit(x) = log(x/(1--x))} for \eqn{x \in (0,1)}{x in (0,1)}.
}

\value{
  Numeric value on requested scale.
}

\note{

\code{logit} is equivalent to \code{\link{qlogis}}, and \code{invlogit} is equivalent to \code{\link{plogis}} (both \R functions in the \pkg{stats} package).
\code{logit} and \code{invlogit} are used in \pkg{secr} because they are slightly more robust to bad input, and their names are more memorable!

}

\examples{
logit(0.5)
invlogit(logit(0.2))
}


