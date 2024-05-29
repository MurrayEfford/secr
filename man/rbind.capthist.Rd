\name{rbind.capthist}
\alias{rbind.capthist}
\alias{MS.capthist}
\alias{append.capthist}

\title{ Combine capthist Objects }

\description{

  Form a single \code{capthist} object from two or more compatible \code{capthist} objects.

}

\usage{

  MS.capthist(...)
  append.capthist(..., synchronous = TRUE)
  
  \method{rbind}{capthist}(..., renumber = TRUE, pool = NULL, verify = TRUE)

}

\arguments{
  \item{\dots}{ one or more \code{capthist} objects or lists of
      \code{capthist} objects }
  \item{synchronous}{logical; if TRUE occasions are assumed to coincide}
  \item{renumber}{ logical, if TRUE assigns new composite individual ID }
  \item{pool}{ list of vectors of session indices or names }
  \item{verify}{ logical, if TRUE the output is checked with \code{\link{verify}} }
}

\details{

  \code{MS.capthist} concatenates the sessions in the input objects as
  one multi-session capthist object. Each session may use a different
  detector array (traps) and a different number of sampling
  occasions. Session names are derived implicitly from the inputs, or
  may be given explicitly (see Examples); if any name is duplicated, all
  will be replaced with sequential integers. The \dots argument may include
  lists of single-session capthist objects.
  
  \code{append.capthist} constructs a single-session capthist assuming the inputs 
  are independent capthist objects. Individuals and detectors are renumbered. Occasions are assumed to be synchronous by default. This is an alternative to rbind when the detector array varies.
  
  The \code{rbind} method for capthist objects is used to pool capture data 
  from more than one
  session into a single session. The number of rows in the output
  session is the sum of the number of rows in the input sessions
  (i.e. each animal appears in only one session). Sessions to be pooled with 
  \code{rbind} must have the same number of capture occasions and use the 
  same detectors (traps). At present there is no function to pool capthist 
  data from different detector arrays. For this it is recommended that you 
  merge the input files and rebuild the capthist object from scratch.

  For \code{rbind.capthist}, the \dots argument may be \enumerate{
      \item A series of single-session capthist objects, which are
    pooled to form one new single-session object, or
      \item One multi-session capthist object, when the components of
    `pool' are used to define combinations of old sessions; e.g. \code{pool =
    list(A=1:3, B=4:5)} produces an object with two sessions (named `A'
    and `B') from 5 old ones. If \code{pool = NULL} (the default) then all the
    sessions are pooled to form one single-session capthist object.
    }

    The names of arguments other than \dots should be given in full. If
\code{renumber = TRUE} (the default), the session name will be prepended
to the animal ID before pooling: animals 1, 2 and 3 in Session A will
become A.1, A.2 and A.3, while those in Session B become B.1, B.2 and
B.3. This ensures that each animal has a unique ID. If \code{renumber =
FALSE}, the animal IDs will not change.

  Other attributes (xy, signal) are handled appropriately. If the signal
    threshold (attribute `cutval') differs among sessions, the maximum is
    used and detections of lower signal strength are discarded.
    
The use of \code{rbind.capthist} to concatenate sessions is now deprecated: 
use \code{MS.capthist}. 

Although \code{MS.capthist} looks like an S3 method, it isn't. The full function name must be used. \code{rbind.capthist} became an S3 method in \pkg{secr} 3.1, so it is called as \code{rbind} alone.

}

\value{

For \code{MS.capthist}, a multi-session object of class `capthist' with
number of sessions equal to the number of sessions in the objects in
\dots.

For \code{append.capthist}, a single-session of class `capthist'. When `synchronous' the dimensions are (\eqn{\sum n_i, \mbox{max}~ S_i, \sum K_i}) where (\eqn{n_i, S_i, K_i}) are the numbers of individuals, occasions and detectors in the i-th input. Otherwise the dimensions are (\eqn{\sum n_i, \sum S_i, \sum K_i}). An index to the original component is saved in the traps covariate `sub'.

For \code{rbind.capthist}, either an object of class `capthist' with one 
session formed by pooling the sessions in the input objects, or a 
capthist object with more than one session, each formed by pooling 
groups of sessions defined by the `pool' argument. Covariate columns that appear in
    all input sessions are retained in the output. 
}

\seealso{ \code{\link{capthist}}, \code{\link{subset.capthist}} }

\examples{

## extend a multi-session object
## we fake the 2010 data by copying from 2005
## note how we name the appended session
fakeCH <- ovenCH[["2005"]]
MS.capthist(ovenCH, "2010" = fakeCH)

## simulate sessions for 2-part mixture
temptrap <- make.grid(nx = 8, ny = 8)
temp1 <- sim.capthist(temptrap,
    detectpar = list(g0 = 0.1, sigma = 40))
temp2 <- sim.capthist(temptrap,
    detectpar = list(g0 = 0.2, sigma = 20))

## concatenate sessions
temp3 <- MS.capthist(large.range = temp1, small.range = temp2)
summary(temp3)
## session-specific movement statistic
RPSV(temp3)

## pool sessions
temp4 <- rbind(temp1, temp2)
summary(temp4)
RPSV(temp4)

## compare mixture to sum of components
## note `detectors visited' is not additive for 'multi' detector
## nor is `detectors used'
(summary(temp1)$counts +  summary(temp2)$counts) -
    summary(temp4)$counts

\dontrun{

## compare two different model fits 
tempfit3 <- secr.fit(temp3, CL = TRUE, buffer = 150, model = list
    (g0 ~ session, sigma ~ session), trace = FALSE)
predict(tempfit3)

## if we can tell which animals had large ranges...
covariates(temp4) <- data.frame(range.size = rep(c("large",
    "small"), c(nrow(temp1), nrow(temp2))))
tempfit4 <- secr.fit(temp4, CL = TRUE, buffer = 150, model = list
    (g0 ~ range.size, sigma ~ range.size), trace = FALSE)
predict(tempfit4, newdata = data.frame(range.size = c("large",
    "small")))

## polygon data
pol1 <- make.poly()
pol2 <- make.poly(x = c(50,50,150,150))
ch1 <- sim.capthist(pol1, popn = list(D = 30), detectfn = 'HHN', 
   detectpar = list(lambda0 = 0.3))
ch2 <- sim.capthist(pol2, popn = list(D = 30), detectfn = 'HHN', 
   detectpar = list(lambda0 = 0.3))
plot(ch1); plot(pol2, add = TRUE); plot(ch2, add = TRUE)

}

}

\keyword{ manip }

