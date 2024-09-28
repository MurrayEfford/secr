\name{plot.secr}
\alias{plot.secr}
\alias{plot.secrlist}
\alias{detectfnplot}
\alias{attenuationplot}

\title{ Plot Detection Functions }

\description{
  Plot detection functions using estimates of parameters in an secr object, or as provided by the user.
}

\usage{
\method{plot}{secr}(x, newdata = NULL, add = FALSE,
    sigmatick = FALSE, rgr = FALSE, limits = FALSE, alpha = 0.05,
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...)

\method{plot}{secrlist}(x, newdata = NULL, add = FALSE,
    sigmatick = FALSE, rgr = FALSE, limits = FALSE, alpha = 0.05,
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...,
    overlay = TRUE)

detectfnplot (detectfn, pars, details = NULL, add = FALSE,
    sigmatick = FALSE, rgr = FALSE, hazard = FALSE, xval = 0:200, ylim = NULL,
    xlab = NULL, ylab = NULL, ...)

attenuationplot (pars, add = FALSE, spherical = TRUE,
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...)

}

\arguments{

  \item{x}{ an \code{secr} object }
  \item{newdata}{ dataframe of data to form estimates }
  \item{add}{ logical to add curve(s) to an existing plot }
  \item{sigmatick}{ logical; if TRUE the scale parameter sigma is shown by a vertical line }
  \item{rgr}{ logical; if TRUE a scaled curve r.g(r) is plotted instead of g(r) }
  \item{hazard}{logical; if TRUE the hazard of detection is plotted instead of probability}
  \item{limits}{ logical; if TRUE pointwise confidence limits are drawn }
  \item{alpha}{ alpha level for confidence intervals }
  \item{xval}{ vector of distances at for which detection to be plotted }
  \item{ylim}{ vector length 2 giving limits of y axis }
  \item{xlab}{ label for x axis }
  \item{ylab}{ label for y axis }
  \item{\dots}{ arguments to pass to \code{lines} }
  \item{overlay}{ logical; if TRUE then automatically \code{add = TRUE}
  for plots after the first }
  \item{detectfn}{integer code or character string for shape of
        detection function 0 = halfnormal etc. -- see \link{detectfn} }
  \item{pars}{ list, vector or matrix of parameter values }
  \item{details}{ list of ancillary parameters }
  \item{spherical}{ logical for whether to include spherical spreading term }

}
\details{

\code{newdata} is usually NULL, in which case one curve is plotted for
each session and group. Otherwise, \code{predict.secr} is used to form
estimates and plot a curve for each row in \code{newdata}.

If axis labels are not provided they default to `Distance (m)' and `Detection probability' or `Detection lambda'.

\code{detectfnplot} is an alternative in which the user nominates the
type of function and provides parameter values. \code{pars} maybe a list
as from \code{\link{detectpar}}; it is first coerced to a numeric vector
with \code{unlist}. Parameter values must be in the expected order
(e.g. g0, sigma, z). If \code{pars} is a matrix then a separate
curve is plotted with the parameter values in each row.

For \code{detectfnplot} the signal threshold parameters `cutval' and
`spherical' should be provided in \code{details} (see examples).

Approximate confidence limits for g(r) are calculated using a numerical
first-order delta-method approximation to the standard error at each
\code{xval}. The distribution of g(r) is assumed to be normal on the logit scale for non-hazard functions (detectfn 0:13). For hazard detection functions (detectfn 14:18) the hazard is assumed (from version 3.1.1) to be distributed normally on the log scale. Limits are back-transformed to the probability scale g(r).

\code{attenuationplot} plots the expected decline in signal strength
with distance, given parameters \eqn{\beta_0}{beta0} and
\eqn{\beta_1}{beta1} for a log-linear model of sound attenuation.

}

\value{

\code{plot.secr} invisibly returns a dataframe of the plotted values (or
a list of dataframes in the case that \code{newdata} has more than one
row).

}

\seealso{ \code{\link{Detection functions}}, \code{\link{plot}}, \code{\link{secr}} }

\examples{

plot (secrdemo.b, xval = 0:100, ylim = c(0, 0.4))
## Add recapture probability
plot (secrdemo.b, newdata = data.frame(b = 1), add = TRUE,
    col = "red")

## signal strength detection: 70dB at source, attenuation
## 0.3dB/m, sdS 5dB; detection threshold 40 dB.
detectfnplot (detectfn = 10, c(70, -0.3, 5), details =
    list(cutval = 40))

## add a function with louder source and spherical spreading...
detectfnplot (detectfn = 11, c(110, -0.3, 5), details =
    list(cutval = 40), add = TRUE, col = "red")

## matching sound attenuation curves; `spherical-only' dashed line
attenuationplot (c(70, -0.3), spherical = FALSE, ylim=c(-10,110))
attenuationplot (c(110, 0), spherical = TRUE, add=TRUE, lty=2)
attenuationplot (c(110, -0.3), spherical = TRUE, add = TRUE,
    col = "red")

}


