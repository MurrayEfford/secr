\name{ms}
\alias{ms}
\alias{ms.default}
\alias{ms.mask}
\alias{ms.secr}

\title{
Multi-session Objects
}

\description{
  Logical function to distinguish objects that span multiple sessions
}

\usage{
  \method{ms}{default}(object, ...)
  \method{ms}{mask}(object, ...)
  \method{ms}{secr}(object, ...)
 }

\arguments{
  \item{object}{any object}
  \item{\dots}{other arguments (not used) }
}
\details{
  
The test applied varies with the type of object. The default method uses
\code{inherits(object, "list")}.

}

\value{
logical, TRUE if \code{object} contains data for multiple sessions
}

\seealso{
\code{\link{capthist}}, \code{\link{mask}}, \code{\link{secr.fit}}
}

\examples{
ms(ovenCH)
ms(ovenbird.model.1)
ms(ovenCH[[1]])
}

