\name{nontarget}
\alias{nontarget}
\alias{interference}

\title{ Non-target Data }

\description{

Non-target detections and interference events may be recorded in a binary detector 
x occasion matrix that is attached as an attribute to a single-session capthist object. 

Matrix entries should be zero for trap x occasion combinations that were not used 
(see \link{usage}) or for which there is a corresponding detection.

Models fitted by \pkg{secr} make no use of these data. They may be used in \pkg{ipsecr}.

From \pkg{secr} 4.5.6, a summary of nontarget data is reported by the 
\code{\link{summary}} method for capthist objects, and the \code{\link{verify}} 
method reports clashes between detections and nontarget data.

}

\examples{

set.seed(123)
ch <- captdata

# traps that caught something
caught <- t(apply(ch, 2:3, sum))

# construct artificial nontarget data
# (positive for half the traps that caught nothing)
attr(ch, 'nontarget') <- (1-caught) * (runif(500)>0.5)

head(caught)
head(attr(ch, 'nontarget'))

# the summary method recognises the 'nontarget' attribute
summary(ch)$nontarget

}

