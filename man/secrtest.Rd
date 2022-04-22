\name{secrtest}
\alias{secrtest}
\alias{print.secrtest}
\alias{plot.secrtest}

\title{ Goodness-of-fit Test Results }

\description{

  S3 class for results from \code{\link{secr.test}}.

}

\usage{

\method{print}{secrtest}(x, terse = TRUE, ...)
\method{plot}{secrtest}(x, stat, ...)

}

\arguments{

  \item{x}{secrtest object from \code{\link{secr.test}} }
  \item{terse}{logical; if TRUE only p values are displayed}
  \item{stat}{character; names of statistics to plot (default: all)}
  \item{...}{other arguments passed to hist by \code{plot.secrtest}}

}

\details{
    
An `secrtest' object is output from \code{\link{secr.test}}.

\code{plot.secrtest} plots a histogram of the simulated values.

If \code{plot.secrtest} is applied to an object with more than one
statistic then multiple plots are produced, so a multi-figure layout
should be prepared (\code{par(mfrow = c(1,2))} for 2 plots side by
side). Include the \code{hist} argument \code{main = ''} to suppress the
ugly plot labels, and ensure each statistic is named by \code{statfn} so
that the x-axis is labelled correctly (See the Examples in help for
\code{\link{secr.test}}).

}

\examples{

\dontrun{

tmp <- secr.test(ovenbird.model.1)
if (inherits(tmp, 'secrtest')) {
    tmp  ## terse print
    print(tmp, terse = FALSE)
    par(mfrow = c(1,5))
    plot(tmp, main = '', xlim=c(0,1), breaks=seq(0,1,0.05))
    par(mfrow = c(1,1))    ## reset to default
}

}

}

\seealso{
  
\code{\link{secr.test}}

}

\keyword{ classes }
