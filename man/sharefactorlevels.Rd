\name{shareFactorLevels}
\alias{shareFactorLevels}

\title{Fix Inconsistent Covariates}

\description{

Factor covariates can give trouble in multi-session models if the levels differ among sessions. A warning is provided by \code{verify.capthist}. This function forces factor covariates to use the same levels. 

}

\usage{

shareFactorLevels(object, columns = NULL, stringsAsFactors = TRUE)

}

\arguments{
  \item{object}{multi-session capthist object or list of traps or mask objects}
  \item{columns}{indices of columns to fix (default all)}
  \item{stringsAsFactors}{logical; if TRUE then character columns are converted to factor}
}

\details{

Factor-valued covariates are coerced to use the same set of levels for each session of a multi-session capthist object or each component of a list of traps or masks. The combined level set is the union of all levels in separate sessions. The order of levels follows the default in \code{\link{factor}} (alphabetical according to current locale).

Setting \code{stringsAsFactors = TRUE} causes character-valued columns to be converted to factors.

}

\value{

An object of the same class as input. A single-session object is passed unchanged except for possible conversion of character values to factor (\code{stringsAsFactors = TRUE}).

}

\seealso{

\code{\link{verify.capthist}}

}

