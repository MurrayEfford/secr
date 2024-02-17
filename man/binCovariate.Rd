\name{binCovariate}
\alias{binCovariate}

\title{Add Binned Covariate}

\description{

Forms a new covariate, replacing values of an old covariate by the central value of equal-width bins.

}

\usage{

binCovariate(object, covname, width)

}

\arguments{

  \item{object}{secr object with covariates attribute (capthist, traps, mask)}
  \item{covname}{character name of covariate}
  \item{width}{numeric bin width}

}

\details{

The name of the new covariate is \code{paste0(covname, width)}.

Fails if covariate not found or is not numeric or there is already a covariate with the new name.

Multi-session objects are handled appropriately. 

}

\value{

Object of the same class as the input with new covariate.

}

\seealso{

\code{\link{covariates}}, \link{skink}

}

\examples{

# bin values of skink snout-vent length (mm)
infraCH <- binCovariate (infraCH, 'SVL', 5)
table(covariates(infraCH[[1]])$SVL5)

# bin values of trap covariate 'HtBrack' (height of bracken, cm)
traps(infraCH) <- binCovariate(traps(infraCH), "HtBrack", 20)
table(covariates(traps(infraCH)[[1]])$HtBrack20)

}
