\name{intervals}
\alias{intervals}
\alias{intervals<-}
\alias{sessionlabels}
\alias{sessionlabels<-}

\title{ Work with Open Population data }

\description{
Functions for data manipulation 
}

\usage{

intervals(object, ...)
intervals(object) <- value
sessionlabels(object, ...)
sessionlabels(object) <- value

}

\arguments{
  \item{object}{ capthist object }
  \item{value}{ vector of intervals or primary session labels}
  \item{\dots}{ other arguments (not used) }
}

\details{

  \code{intervals} extracts the `interval' attribute if it exists.

  The attribute `intervals' is set automatically by the secr function
  \code{\link{join}}.
  
 \code{sessionlabels} provides session names for the primary sessions encoded 
 in a ``single-session'' capthist object (e.g., the result of \code{join}) 
 that has an intervals attribute. The names are used by some summary functions 
 in the package \pkg{openCR} (M. Efford unpubl.) (\code{m.array}, \code{JS.counts}).

The function \code{\link{session}} has a different purpose: labelling
sessions in a multi-session capthist object. However, \code{session}
names of multi-session input are used automatically by \code{\link{join}} to construct the
\code{sessionlabels} attribute of the resulting single-session object.

}

\value{
  For \code{intervals}, a numeric vector of time intervals, one less than the number of occasions (secondary sessions). 
  
  For \code{sessionlabels}, a character vector of primary session names.

}

\note{

  There is a naming conflict with the intervals function in \pkg{nlme}.

}

\examples{

singlesessionCH <- join(ovenCH)
intervals(singlesessionCH)
sessionlabels(singlesessionCH)

}

