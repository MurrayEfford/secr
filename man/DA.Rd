\name{BUGS}
\alias{read.DA}
\alias{write.DA}

\title{
Convert Data To Or From BUGS Format
}

\description{
Convert data between `capthist' and BUGS input format.
}
\usage{

read.DA(DAlist, detector = "polygonX", units = 1, session = 1,
    Y = "Y", xcoord = "U1", ycoord = "U2", xmin = "Xl",
    xmax = "Xu", ymin = "Yl", ymax = "Yu", buffer = "delta",
    verify = TRUE)

write.DA(capthist, buffer, nzeros = 200, units = 1)

}

\arguments{
  \item{DAlist}{ list containing data in BUGS format }
  \item{detector}{ character value for detector type: `polygon' or `polygonX' }
  \item{units}{ numeric for scaling output coordinates}
  \item{session}{numeric or character label used in output}
  \item{Y}{character, name of binary detection history matrix (animals x occasions)}
  \item{xcoord}{character, name of matrix of x-coordinates for each detection in \code{Y}}
  \item{ycoord}{character, name of matrix of y-coordinates for each detection in \code{Y}}
  \item{xmin}{ character, name of coordinate of state space boundary}
  \item{xmax}{ character, name of coordinate of state space boundary}
  \item{ymin}{ character, name of coordinate of state space boundary}
  \item{ymax}{ character, name of coordinate of state space boundary}
  \item{buffer}{ see Details}
  \item{verify}{ logical if TRUE then the resulting capthist object is
    checked with \code{\link{verify}}}
  \item{capthist}{\code{\link{capthist}} object }
  \item{nzeros}{level of data augmentation (all-zero detection histories)}
}

\details{

  Data for OpenBUGS or WinBUGS called from R using the package
  \pkg{R2WinBUGS} (Sturtz et al. 2005) take the form of an \R list.

  These functions are limited at present to binary data from a square
  quadrat such as used by Royle and Young (2008). Marques et al. (2011)
  provide an \R function \code{create.data()} for generating simulated
  datasets of this sort (see \code{\link{sim.capthist}} for equivalent
  functionality).

  When reading BUGS data --
  
  The character values \code{Y}, \code{xcoord}, \code{ycoord},
  \code{xmin} etc. are used to locate the data within \code{DAlist},
  allowing for variation in the input names.

  The number of sampling occasions is taken from the number of columns
  in \code{Y}. Each value in \code{Y} should be 0 or 1. Coordinates may
  be missing
  
  A numeric value for \code{buffer} is the distance (in the original
  units) by which the limits Xl, Xu etc. should be shrunk to give the
  actual plot limits.  If \code{buffer} is character then a component of
  \code{DAlist} contains the required numeric value.

  Coordinates in the output will be \emph{multiplied by} the scalar
  \code{units}.

  Augmentation rows corresponding to `all-zero' detection histories in
  \code{Y}, \code{xcoord}, and \code{ycoord} are discarded.

  When writing BUGS data --

  Null (all-zero) detection histories are added to the matrix of
  detection histories \code{Y}, and missing (NA) rows are added to the
  coordinate matrices \code{xcoord} and \code{ycoord}.

  Coordinates in the output will be \emph{divided by} the scalar
  \code{units}.

}

\value{

  For \code{read.DA}, an object of class `capthist'.
  
  For \code{write.DA}, a list with the components

\tabular{lll}{
Xl \tab left edge of state space \cr
Xu \tab right edge of state space \cr
Yl \tab bottom edge of state space \cr
Yu \tab top edge of state space \cr
delta \tab buffer between edge of state space and quadrat \cr
nind \tab number of animals observed \cr
nzeros \tab number of added all-zero detection histories \cr
T \tab number of sampling occasions \cr
Y \tab binary matrix of detection histories (dim = c(nind+nzeros, T)) \cr
U1 \tab matrix of x-coordinates, dimensioned as Y \cr
U2 \tab matrix of y-coordinates, dimensioned as Y \cr
}

U1 and U2 are `NA' where animal was not detected.

}

\references{

Marques, T. A., Thomas, L. and Royle, J. A. (2011) A hierarchical model
for spatial capture--recapture data: Comment. \emph{Ecology} \bold{92},
526--528.

Royle, J. A. and Young, K. V. (2008) A hierarchical model for spatial
capture--recapture data. \emph{Ecology} \bold{89}, 2281--2289.

Sturtz, S., Ligges, U. and Gelman, A. (2005) R2WinBUGS: a package for
    running WinBUGS from R. \emph{Journal of Statistical Software}
    \bold{12}, 1--16.
}

\seealso{

  \code{\link{hornedlizardCH}}, \code{\link{verify}}, \code{\link{capthist}}

}

\examples{

write.DA (hornedlizardCH, buffer = 100, units = 100)

## In this example, the input uses Xl, Xu etc.
## for the limits of the plot itself, so buffer = 0.
## Input is in hundreds of metres.
## First, obtain the list lzdata
olddir <- setwd (system.file("extdata", package="secr"))
source ("lizarddata.R")
setwd(olddir)
str(lzdata)
## Now convert to capthist
tempcapt <- read.DA(lzdata, Y = "H", xcoord = "X",
    ycoord = "Y", buffer = 0, units = 100)

\dontrun{

plot(tempcapt)
secr.fit(tempcapt, trace = FALSE)
## etc.

}

}

