\name{detector}
\alias{detector}
\alias{detector<-}
\alias{polygon}
\alias{polygonX}
\alias{transect}
\alias{transectX}
\alias{single}
\alias{multi}
\alias{proximity}
\alias{count}
\alias{capped}

% \alias{signal}

\title{ Detector Type }

\description{
  Extract or replace the detector type. 
}

\usage{
detector(object, \dots)
detector(object) <- value
}

\arguments{
  \item{object}{ object with `detector' attribute e.g. \code{traps} }
  \item{value}{ character string for detector type }
  \item{\dots}{ other arguments (not used) }
}

\details{

Valid detector types are `single', `multi', `proximity', `count', `capped',
`signal', `polygon', `transect', `polygonX', and `transectX'. The
detector type is stored as an attribute of a \code{traps} object.
Detector types are mostly described by Efford et al. (2009a,b; see also
\href{../doc/secr-overview.pdf}{secr-overview.pdf}). Polygon and transect detector types are
for area and linear searches as described in
\href{https://www.otago.ac.nz/density/pdfs/secr-polygondetectors.pdf}{secr-polygondetectors.pdf} and Efford (2011). The `signal'
detector type is used for acoustic data as described in
\href{https://www.otago.ac.nz/density/pdfs/secr-sound.pdf}{secr-sound.pdf}.

The `capped' detector type refers to binary proximity data in which no more than one individual may be detected at a detector on any occasion. The type is partially implemented in \pkg{secr} 3.1.1: data may be simulated and manipulated, but for model fitting these are treated as proximity data by \code{secr.fit()}.
}

\value{
  character string for detector type
}

\seealso{ 
  \code{\link{traps}}, \code{\link{RShowDoc}}
}

\references{

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area searches. \emph{Ecology}
\bold{92}, 2202--2207.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009a) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009b) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\examples{
## Default detector type is "multi"
temptrap <- make.grid(nx = 6, ny = 8)
detector(temptrap) <- "proximity"
summary(temptrap)
}

\keyword{ models }

