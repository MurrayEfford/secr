\name{plot.MCgof}
\alias{plot.MCgof}
\alias{print.MCgof}
\alias{summary.MCgof}

\title{Methods for MCgof Objects}

\description{

Plot, summary and print methods for MCgof objects.

}

\usage{

\method{plot}{MCgof}(x, counts = 'all', overlay = NULL, maxT = NULL, 
    main = NULL, cex = 0.9, ...)

\method{summary}{MCgof}(object, ...) 

\method{print}{MCgof}(x, ...) 

}

\arguments{

  \item{x}{MCgof object}
  \item{counts}{character vector of marginal counts for which statistics are to be plotted}
  \item{overlay}{MCgof object}
  \item{maxT}{numeric maximum plotted value of statistic}
  \item{main}{character vector of labels (see Details)}
  \item{cex}{numeric size of labels and points}
  \item{\dots}{other arguments passed by the plot method to \code{\link{points}} for plotting overlay}
  \item{object}{MCgof object}

}

\details{

We start with a 3-D capthist array with dimensions corresponding to individuals (i), occasions (j) and detectors (k). The possible marginal counts for the default `statfn' in \code{\link{MCgof}} are designated --

\tabular{lll}{
Count \tab Margin \tab Cell value \cr
yik \tab individual x detector \tab \eqn{y_{ik} = \sum_j y_{ijk}} \cr
yi \tab individual \tab \eqn{y_i = \sum_j \sum_k y_{ijk}} \cr
yk \tab detector \tab \eqn{y_k = \sum_j \sum_i y_{ijk}} \cr
}

The plot method displays a scatterplot of discrepancies for observed and simulated data (one point per replicate) (Gelman et al. 1996). 

If `overlay' is provided then the results are overlaid on the initial plot. Points should be distinguished by specifying a different colour (col) or symbol (pch) with the \dots argument.

`main' is a vector of labels used as headers; the names should include all components of `statfn'. Setting main = "" suppresses headers.

}

\value{

The summary method returns a matrix of values in which the columns correspond to the different statistics (default yik, yi, yk) and the rows are

\itemize{
\item median discrepancy Tobs
\item median discrepancy Tsim
\item proportion Tobs>Tsim
\item number of valid results
}

}

\seealso{

\code{\link{MCgof}}

}

\references{

Gelman, A., Meng, X.-L., and Stern, H. (1996) Posterior predictive assessment of 
model fitness via realized discrepancies. \emph{Statistica Sinica} \bold{6}, 733--807.

}
