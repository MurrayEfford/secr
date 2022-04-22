\name{reduce.capthist}
\alias{reduce.traps}
\alias{reduce.capthist}

\title{ Combine Occasions Or Detectors }

\description{

Use these methods to combine data from multiple occasions or multiple detectors in a
\code{capthist} or \code{traps} object, creating a new data set and possibly converting
between detector types.

}

\usage{

\method{reduce}{traps}(object, newtraps = NULL, newoccasions = NULL, span = NULL, 
    rename = FALSE, newxy = c('mean', 'first'), \dots)

\method{reduce}{capthist}(object, newtraps = NULL, span = NULL, rename =
    FALSE, newoccasions = NULL, by = 1, outputdetector = NULL,
    select = c("last","first","random"), dropunused = TRUE, verify
    = TRUE, sessions = NULL, \dots)
}

\arguments{

  \item{object}{ \code{traps or capthist} object }

  \item{newtraps}{ list in which each component is a vector of
  subscripts for detectors to be pooled }

  \item{newoccasions}{ list in which each component is a vector of
  subscripts for occasions to be pooled }

  \item{span}{ numeric maximum span in metres of new detector }

  \item{rename}{logical; if TRUE the new detectors will be numbered from
1, otherwise a name will be constructed from the old detector names}

   \item{newxy}{ character; coordinates when detectors grouped with `newtraps'}  

  \item{by}{ number of old occasions in each new occasion }

  \item{outputdetector}{character value giving \link{detector} type for output (defaults to input) }

  \item{select}{ character value for method to resolve conflicts }

  \item{dropunused}{ logical, if TRUE any never-used detectors are dropped }

  \item{verify}{ logical, if TRUE the \code{verify} function is applied to the output }

  \item{sessions}{ vector of session indices or names (optional) }

  \item{\dots}{ other arguments passed by reduce.capthist to reduce.traps, 
  or by reduce.traps to hclust }

}

\details{

 
  reduce.traps --
  
Grouping may be specified explicitly via \code{newtraps}, or
implicitly by \code{span}. 

If \code{span} is specified a clustering of detector sites will be
performed with \code{\link{hclust}} and detectors will be assigned to
groups with \code{\link{cutree}}. The default algorithm in \code{hclust}
is complete linkage, which tends to yield compact, circular clusters;
each will have diameter less than or equal to \code{span}.

\code{newxy = 'first'} selects the coordinates of the first detector in a 
group defined by `newtraps', rather then the average of all detectors in group.

reduce.capthist --
  
The first component of \code{newoccasions} defines the columns of
\code{object} for new occasion 1, the second for new occasion 2, etc. If
\code{newoccasions} is NULL then all occasions are output. Subscripts in a
component of \code{newoccasions} that do not match an occasion in the input
are ignored. When the output detector is one of the trap types
(`single', `multi'), reducing capture occasions can result in locational
ambiguity for individuals caught on more than one occasion, and for
single-catch traps there may also be conflicts between individuals at
the same trap. The method for resolving conflicts among `multi'
detectors is determined by \code{select} which should be one of `first',
`last' or `random'. With `single' detectors \code{select} is ignored and
the method is: first, randomly select* one trap per animal per day;
second, randomly select* one animal per trap per day; third, when
collapsing multiple days use the first capture, if any, in each
trap.

Usage data in the \code{traps} attribute are also pooled if present;
usage is summed over contributing occasions and detectors. If there is
no 'usage' attribute in the input, and \code{outputdetector} is one of
'count', 'polygon', 'transect' and 'telemetry', a homogeneous (all-1's)
'usage' attribute is first generated for the input.

* i.e., in the case of a single capture, use that capture; in the case of
multiple `competing' captures draw one at random.

If \code{newoccasions} is not provided then old occasions are grouped into
new occasions as indicated by the \code{by} argument. For example, if
there are 15 old occasions and \code{by = 5} then new occasions will be
formed from occasions 1:5, 6:10, and 11:15. A warning is given when the
number of old occasions is not a multiple of \code{by} as then the final
new occasion will comprise fewer old occasions.

\code{dropunused = TRUE} has the possibly unintended effect of dropping
whole occasions on which there were no detections.

A special use of the \code{by} argument is to combine all occasions into
one for each session in a multi-session dataset. This is done by setting
\code{by = "all"}.

\code{reduce.capthist} may be used with non-spatial capthist objects 
(NULL 'traps' attribute) by setting \code{verify = FALSE}.

}

\value{

  reduce.traps --

  An object of class traps with detectors combined according to
  \code{newtraps} or \code{span}. The new object has an attribute
  `newtrap', a vector of length equal to the original number of
  detectors. Each element in newtrap is the index of the new detector to
  which the old detector was assigned (see Examples).

  The object has no clusterID or clustertrap attribute.
  
  reduce.capthist --
  
An object of class capthist with number of occasions (columns) equal to
\code{length(newoccasions)}; detectors may simulataneously be aggregated
as with \code{reduce.traps}. The detector type is inherited from \code{object}
unless a new type is specified with the argument
\code{outputdetector}. 

}

\seealso{
  
  \code{\link{capthist}}, \code{\link{subset.capthist}},
  \code{\link{discretize}}, \code{\link{hclust}}, \code{\link{cutree}}
  
}

\section{Warning}{The argument named `columns' was renamed to
  `newoccasions' in version 2.5.0, and arguments were added to
  reduce.capthist for the pooling of detectors. Old code should work as
  before if all arguments are named and `columns' is changed. }

\note{

  The reduce method may be used to re-assign the detector type (and
  hence data format) of a capthist object without combining occasions or
  detectors. Set the \code{object} and \code{outputdetector} arguments
  and leave others at their default values.
  
  Automated clustering can produce unexpected outcomes. In particular,
  there is no guarantee that clusters will be equal in size. You should
  inspect the results of reduce.traps especially when using \code{span}.

  \code{reduce.traps} is not implemented for polygons or transects.

  The function \code{\link{discretize}} converts polygon data to
  point-detector (multi, proximity or count) data.
  
}

\examples{
tempcapt <- sim.capthist (make.grid(nx = 6, ny = 6), nocc = 6)
class(tempcapt)

pooled.tempcapt <- reduce(tempcapt, newocc = list(1,2:3,4:6))
summary (pooled.tempcapt)

pooled.tempcapt2 <- reduce(tempcapt, by = 2)
summary (pooled.tempcapt2)

## collapse multi-session dataset to single-session 'open population'
onesess <- join(reduce(ovenCH, by = "all"))
summary(onesess)

# group detectors within 60 metres
plot (traps(captdata))
plot (reduce(captdata, span = 60), add = TRUE)

# plot linking old and new
old <- traps(captdata)
new <- reduce(old, span = 60)
newtrap <- attr(new, "newtrap")
plot(old, border = 10)
plot(new, add = TRUE, detpar = list(pch = 16), label = TRUE)
segments (new$x[newtrap], new$y[newtrap], old$x, old$y)

\dontrun{

# compare binary proximity with collapsed binomial count
# expect TRUE for each year
for (y in 1:5) {
    CHA <- abs(ovenCHp[[y]])   ## abs() to ignore one death
    usage(traps(CHA)) <- matrix(1, 44, ncol(CHA))
    CHB <- reduce(CHA, by = 'all', output = 'count')
    # summary(CHA, terse = TRUE)
    # summary(CHB, terse = TRUE)
    fitA <- secr.fit(CHA, buffer = 300, trace = FALSE)
    fitB <- secr.fit(CHB, buffer = 300, trace = FALSE, binomN = 1, biasLimit = NA)
    A <- predict(fitA)[,-1] 
    B <- predict(fitB)[,-1]
    cat(y, ' ', all(abs(A-B)/A < 1e-5), '\n')
}
## multi-session fit
## expect TRUE overall
CHa <- ovenCHp
for (y in 1:5) {
    usage(traps(CHa[[y]])) <- matrix(1, 44, ncol(CHa[[y]]))
    CHa[[y]][,,] <- abs(CHa[[y]][,,])
}
CHb <- reduce(CHa, by = 'all', output = 'count')
summary(CHa, terse = TRUE)
summary(CHb, terse = TRUE)
fita <- secr.fit(CHa, buffer = 300, trace = FALSE)
fitb <- secr.fit(CHb, buffer = 300, trace = FALSE, binomN = 1, biasLimit = NA)
A <- predict(fita)[[1]][,-1] 
B <- predict(fitb)[[1]][,-1]
all(abs(A-B)/A < 1e-5)

}

}

\keyword{ manip }
