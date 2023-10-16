\name{kfn}
\alias{kfn}

\title{Overlap Index}

\description{

Computes the overlap index of Efford et al. (2016) from various inputs, 
including fitted models.

}

\usage{

kfn(object)

}

\arguments{

  \item{object}{fitted secr model, numeric vector, matrix, dataframe}

}

\details{

\code{kfn} simply computes \eqn{k = \sigma \sqrt D / 100}, where \eqn{\sigma} is the sigma parameter of a fitted halfnormal detection function and \eqn{D} is the corresponding density estimate. The factor of 1/100 adjusts for the units used in \pkg{secr} (sigma in metres; D in animals per hectare).

Input may be in any of these forms

\enumerate{

\item vector with D and sigma in the first and second positions.

\item matrix with each row as in (1)

\item dataframe such as produced by \code{\link{predict.secr}} with rows `D' and `sigma', and column `estimate'.

\item fitted \pkg{secr} model

\item a list of any of the above

}

}

\value{

Numeric vector with elements `D', `sigma' and `k', or a matrix with these columns.

}

\note{

The index should not be taken too literally as a measure of overlap: it represents the overlap expected \emph{if} activity centres are randomly distributed and \emph{if} home ranges have bivariate normal utilisation. Thus it does \emph{not} measure overlap due to social behaviour etc. except as that affects home range size.

The index may be estimated directly using the sigmak parameterization 
(i.e., when sigmak appears in the model for \code{\link{secr.fit}}). 
This provides SE and confidence limits for sigmak (= \eqn{k}). However, 
the directly estimated value of sigmak lacks the unit correction and is 
therefore 100 \eqn{\times} the value from \code{kfn}.

}


\references{

Efford, M. G., Dawson, D. K., Jhala, Y. V. and Qureshi, Q. (2016) 
Density-dependent home-range size revealed by spatially explicit 
capture--recapture. \emph{Ecography} \bold{39}, 676--688.

}

\seealso{

\code{\link{predict.secr}}, \code{\link{secr.fit}}, \code{\link{details}}

}

\examples{

kfn(secrdemo.0)

## compare
\donttest{
fitk <- secr.fit(captdata, model = sigmak~1, buffer = 100, trace = FALSE)
predict(fitk)
}

}
