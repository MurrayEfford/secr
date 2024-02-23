\name{snip}
\alias{snip}

\title{Slice Transect Into Shorter Sections}

\description{

  This function splits the transects in a `transect' or `transectX'
  traps object into multiple shorter sections. The function may also be
  applied directly to a capthist object based on transect data. This
  makes it easy to convert detection data collected along linear
  transects to point detection data (see Example).
  
}

\usage{

snip(object, from = 0, by = 1000, length.out = NULL, keep.incomplete = TRUE, tol = 0.01)

}

\arguments{
  \item{object}{secr `traps' or `capthist' object based on transects}
  \item{from}{numeric starting posiiton (m)}
  \item{by}{numeric length of new transects (m) }
  \item{length.out}{numeric number of new transects, as alternative to `by'}
  \item{keep.incomplete}{logical; if TRUE then initial or terminal
  sections of each original transect that are less than `by' will be
  retained in the output}
  \item{tol}{numeric tolerance for xyontransect (capthist only)}
}

\details{

If a positive \code{length.out} is specified, \code{by} will be computed as
\code{(transectlength(object) - from) / length.out}.

}

\value{

  A `traps' or `capthist' object, according to the input.
  If \code{keep.incomplete == FALSE} animals and detections from the

}

\seealso{
\code{\link{transectlength}}, \code{\link{discretize}}
}

\section{Warning}{

\code{snip} does not work for mark--resight data.

}

\examples{

x <- seq(0, 4*pi, length = 41)
temptrans <- make.transect(x = x*100, y = sin(x)*300)
plot (snip(temptrans, by = 200), markvertices = 1)

\dontrun{

## simulate some captures
tempcapt <- sim.capthist(temptrans, popn = list(D = 2,
   buffer = 300), detectfn = 'HHN', binomN = 0,
   detectpar = list(lambda0 = 0.5, sigma = 50))

## snip capture histories
tempCH <- snip(tempcapt, by = 20)

## collapse from 'transect' to 'count', discarding location within transects
tempCH <- reduce(tempCH, outputdetector = "count")

## fit secr model and examine H-T estimates of density
## fails with detectfn = 'HN'
fit <- secr.fit(tempCH, buffer = 300, CL = TRUE, detectfn = 'HHN', trace = FALSE)
derived(fit)

## also, may split an existing transect into equal lengths
## same result:
plot(snip(temptrans, by = transectlength(temptrans)/10),
    markvertices = 1)
plot(snip(temptrans, length.out = 10), markvertices = 1)

}

}

\keyword{ manip }
