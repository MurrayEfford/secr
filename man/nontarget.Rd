\name{nontarget}
\alias{nontarget}
\alias{nontarget<-}
\alias{interference}

\title{ Non-target Data }

\description{

Non-target detections and interference events may be recorded in a binary detector 
x occasion matrix attached as an attribute to a single-session capthist object, 
or to each component of a multi-session capthist object. 

Models fitted by \pkg{secr} make no use of these data. They may be used in \pkg{ipsecr}.

From \pkg{secr} 4.5.6, a summary of nontarget data is reported by the 
\code{\link{summary}} method for capthist objects, and the \code{\link{verify}} 
method reports clashes between detections and nontarget data.

Extraction and replacement funcions are provided from \pkg{secr} 4.5.7 on.

}

\usage{

nontarget(object, \dots)
nontarget(object) <- value

}

\arguments{

  \item{object}{capthist object}
  \item{\dots}{other argments (not used)}
  \item{value}{ numeric binary matrix (rows = detectors, columns = occasions) }

}

\details{

The order of rows should match the order of detectors in \code{traps(object)}. 
Matrix entries should be zero for trap x occasion combinations that were not used 
(see \link{usage}) or for which there is a corresponding detection.

\code{value} is coerced to a matrix before assignment.

}

\value{

For \code{nontarget()}, a matrix or list of matrices.

}

\seealso{

\code{\link{plot.capthist}}, \code{\link{summary.capthist}}

}


\examples{

set.seed(123)
ch <- captdata

# traps that caught something
caught <- t(apply(ch, 2:3, sum))

# construct artificial nontarget data
# (positive for half the traps that caught nothing)
nontarget(ch) <- (1-caught) * (runif(500)>0.5)

head(caught)
head(nontarget(ch))

# the summary method recognises the 'nontarget' attribute
summary(ch)$nontarget

}

