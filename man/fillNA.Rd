\name{fillNA}
\alias{fillNA}

\title{Fill Missing Covariate Values}

\description{

Missing covariate values are not allowed in a habitat mask. This function replaces 
each missing value with a given value or the closest non-missing value. 
This is not recommended if more than a few cells have missing values. 

}

\usage{

fillNA(mask, cov = NULL, method = c("nearest", "value"), value = NULL, verbose = TRUE)

}

\arguments{
  \item{mask}{habitat mask or similar}
  \item{cov}{character vector of covariate names}
  \item{method}{character}
  \item{value}{character, numeric, or factor value to replace all NA}
  \item{verbose}{logical; if TRUE the replacements are reported}
}

\details{

If there is more than one covariate then each is handled separately. 
If no covariate is named then all covariates are filled.

When method = "value" all NA values are replaced with the given value.

The method argument switches automatically to `value' if value is specified.

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

msk3 <- fillNA(msk, value = factor(2, levels=1:2))

# a call to fillNA shows whether any covariate has missing values
# assign result to tmp to avoid unwanted output
tmp <- fillNA(msk2)

}