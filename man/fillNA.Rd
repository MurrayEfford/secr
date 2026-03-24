\name{fillNA}
\alias{fillNA}

\title{Fill Missing Covariate Values}

\description{

Missing covariate values are not allowed in a habitat mask. This function replaces 
each missing value with the closest non-missing value. It is intended for when a 
few cells have missing values.

}

\usage{

fillNA(mask, cov = NULL, method = "nearest", verbose = TRUE)

}

\arguments{
  \item{mask}{habitat mask or similar}
  \item{cov}{character vector of covariate names}
  \item{method}{character (only 'nearest' implemented)}
  \item{verbose}{logical; if TRUE the replacements are reported}
}

}
\details{

If there is more than one covariate then each is handled separately. 
If no covariate is named then all covariates are filled.

The `method' argument is included to allow for other fill algorithms in the future.

}

\value{

An object with the same class as the input and the requested covariate(s) replaced with non-missing values.}

\seealso{

\link{skink},
\code{\link{addCovariates}}

}

\examples{

# artificial example
msk <- make.mask(traps(captdata), spacing = 20)
covariates(msk) <- data.frame(rand = sample(factor(c(NA,1,2)), nrow(msk), 
    replace = TRUE, prob = c(0.2,0.4,0.4)))

msk2 <- fillNA(msk)

par(mfrow = c(1,2))
plot(msk, cov='rand', dots = FALSE)
plot(msk2, cov='rand', dots = FALSE)

}