\name{trim}
\alias{trim}
\alias{trim.default}
\alias{trim.secr}
\alias{trim.secrlist}

\title{
Drop Unwanted List Components
}

\description{

Drop unwanted components from a \code{list} object, usually to save space.

}

\usage{

\method{trim}{default}(object, drop, keep)
\method{trim}{secr}(object, drop = c("call", "mask", "designD", "designNE", 
    "design", "design0"), keep = NULL)
\method{trim}{secrlist}(object, drop = c("call", "mask", "designD", "designNE", 
    "design", "design0"), keep = NULL)

}

\arguments{ 
    \item{object}{ a list object }
    \item{drop}{ vector identifying components to be dropped } 
    \item{keep}{ vector identifying components to be kept } 
} 

\details{ 

\code{drop} may be a character vector of names or a numeric vector of
indices. If both \code{drop} and \code{keep} are given then the action
is conservative, dropping only components in \code{drop} and not in
\code{keep}.

Be warned that some further operations on fitted secr objects become
impossible once you have discarded the default components.

}

\value{
a list retaining selected components.
}

\examples{

names(secrdemo.0)
names(trim(secrdemo.0))

object.size(secrdemo.0)
object.size(trim(secrdemo.0))

object.size(trim(secrlist(secrdemo.0, secrdemo.b)))

}

\keyword{ manip }

