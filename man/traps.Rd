\name{traps}
\alias{traps object}
\alias{traps}
\alias{traps<-}

\title{ Detector Array }

\usage{ 
traps(object, \dots)
traps(object) <- value
}

\arguments{
\item{object}{ a \code{capthist} object. }
\item{value}{ \code{traps} object to replace previous. }
\item{\dots}{ other arguments (not used). }
}

\description{

An object of class \code{traps} encapsulates a set of detector (trap)
locations and related data. A method of the same name extracts or
replaces the \code{traps} attribute of a \code{capthist} object.

}


\details{

An object of class \code{traps} holds detector (trap) locations as a
data frame of x-y coordinates. Trap identifiers are used as row names.
The required attribute `detector' records the type of detector
("single", "multi" or "proximity" etc.; see \code{\link{detector}} for
more).

Other possible attributes of a \code{traps} object are:

\tabular{ll}{
 \code{\link{spacing}} \tab mean distance to nearest detector \cr
 \code{spacex} \tab \cr
 \code{spacey} \tab  \cr
 \code{\link{covariates}} \tab dataframe of trap-specific covariates \cr
 \code{\link{clusterID}} \tab identifier of the cluster to which each detector belongs  \cr
 \code{\link{clustertrap}} \tab sequence number of each trap within its cluster  \cr  
 \code{\link{usage}} \tab a traps x occasions matrix of effort (may be binary 
 0/1)\cr
 \code{\link{markocc}} \tab integer vector distinguishing marking occasions (1) from 
  sighting occasions (0) \cr
  \code{newtrap} \tab vector recording aggregation of detectors by
 \code{\link{reduce.traps}}\cr

   }

If usage is specified, at least one detector must be `used' (usage non-zero)
    on each occasion.
    
Various array geometries may be constructed with functions such as
\code{\link{make.grid}} and \code{\link{make.circle}}, and these may be
combined or placed randomly with \code{\link{trap.builder}}. 

}

\references{ 

Efford, M. G. (2012) \emph{DENSITY 5.0: software for spatially explicit
capture--recapture}. Department of Mathematics and Statistics,
University of Otago, Dunedin, New Zealand. \url{https://www.otago.ac.nz/density/}.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

}

\note{
  
Generic methods are provided to select rows
(\code{\link{subset.traps}}), combine two or more arrays
(\code{\link{rbind.traps}}), aggregate detectors
(\code{\link{reduce.traps}}), shift an array
(\code{\link{shift.traps}}), or rotate an array
(\code{\link{rotate.traps}}).

The attributes \code{usage} and \code{covariates} may be extracted or
replaced using generic methods of the same name.

}

\seealso{ 

  \code{\link{make.grid}}, \code{\link{read.traps}},
  \code{\link{rbind.traps}}, \code{\link{reduce.traps}},
  \code{\link{plot.traps}}, \code{\link{secr.fit}},
  \code{\link{spacing}}, \code{\link{detector}},
  \code{\link{covariates}}, \code{\link{trap.builder}},
  \code{\link{as.mask}}
  
  }

\examples{

demotraps <- make.grid(nx = 8, ny = 6, spacing = 30)
demotraps    ## uses print method for traps
summary (demotraps)

plot (demotraps, border = 50, label = TRUE, offset = 8, 
    gridlines=FALSE)  

## generate an arbitrary covariate `randcov'
covariates (demotraps) <- data.frame(randcov = rnorm(48))

## overplot detectors that have high covariate values
temptr <- subset(demotraps, covariates(demotraps)$randcov > 0.5)
plot (temptr, add = TRUE, 
    detpar = list (pch = 16, col = "green", cex = 2))  
}


