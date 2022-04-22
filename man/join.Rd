\name{join}
\alias{join}
\alias{unjoin}

\title{ Combine or Split Sessions of capthist Object }

\description{ Make a single-session capthist object from a list of
single-session objects, or a multi-session capthist object.  }

\usage{ 

join(object, remove.dupl.sites = TRUE, tol = 0.001, sites.by.name = FALSE,
    drop.sites = FALSE, intervals = NULL, sessionlabels = NULL, 
    timevaryingcov = NULL) 
    
unjoin(object, intervals, ...)

}

\arguments{
  
  \item{object}{list of single-session objects, or a multi-session
capthist object [\code{join}], or a single-session capthist object [\code{unjoin}]}

  \item{remove.dupl.sites}{ logical; if TRUE then a single record is
retained for each trap site used in multiple input sessions }

  \item{tol}{ absolute distance in metres within which sites are
considered identical }

  \item{sites.by.name}{ logical; if TRUE and \code{remove.dupl.sites} then duplicate sites are inferred from row names rather than x-y coordinates}
  
  \item{drop.sites}{logical; if TRUE then site information is discarded}

  \item{intervals}{ vector of times between sessions (join) or occasions (unjoin; zero indicates same session) }
  
  \item{sessionlabels}{ vector of session names }
  
  \item{timevaryingcov}{ character vector of covariate names }

  \item{...}{ other arguments passed to \code{\link{subset.capthist}} }
}

\details{

  \describe{
\item{\code{join}}{

  The input sessions are assumed to be of the same detector type and to
  have the same attributes (e.g., covariates should be present for all
  or none).
  
  The number of occasions (columns) in the output is equal to the sum of
the number of occasions in each input.

Duplicates may be defined either as sites within a given distance (\code{tol}) or sites with the same name (\code{sites.by.name = TRUE}). Using site names is faster.

For non-spatial analyses it is efficient to drop the third dimension and discard the traps attribute (\code{drop.sites = TRUE}).

A new dataframe of individual covariates is formed using the covariates
for the first occurrence of each animal. 

If \code{timevaryingcov} is given then for each name a new covariate is generated for each session and populated with values observed in that session, or NA if the animal was not detected. A `timevaryingcov' (list) attribute is created that associates each set of new session-specific columns with the corresponding old name, so that it may be used in formulae (see \code{\link{timevaryingcov}}).

Attributes xy and signal are handled appropriately, as is trap usage.
}
\item{\code{unjoin}}{

The input grouping of occasions (columns) into sessions is specified via
\code{intervals}. This is a vector of length one less than the number of
occasions (columns) in \code{object}. Elements greater than zero
indicate a new session.

The \code{intervals} argument may be omitted if \code{object} has a
valid `intervals' attribute, as in the output from \code{join}.
}
}
}
\value{

For \code{join}, a single-session capthist object. The vector attribute `intervals' records the
distinction between occasions that are adjacent in the input (interval =
0) and those that are in consecutive sessions (e.g., interval = 1); `intervals'
has length one less than the number of occasions.

For \code{unjoin}, a multi-session capthist object. Sessions are named
with integers.

}

\note{

Do not confuse \code{unjoin} with \code{\link{split.capthist}} which
splits by row (animal) rather than by column (occasion).

Occasions survive intact; to pool occasions use
\code{\link{reduce.capthist}}.

\code{join} was modified in version 2.9.5 to check whether the
components of `object' all used the same detectors (`traps') (putting
aside differences in usage). If the traps are identical and
remove.dupl.sites = TRUE then the resulting `capthist' uses the common list of
detectors, with a usage attribute formed by concatenating the usage
columns of the input. This is faster than the previous filtering
algorithm using `tol'; the older algorithm is still used if the traps differ. 

Problems may be encountered with large datasets. These may be alleviated by setting sites.by.name = TRUE (if matching sites have matching names, avoiding the need for coordinate matching) or drop.sites = TRUE (if only non-spatial data are required for openCR).

}

\seealso{
  
  \code{\link{MS.capthist}}, \code{\link{rbind.capthist}}
  
}

\examples{

joined.ovenCH <- join (ovenCH)
summary(joined.ovenCH)
attr(joined.ovenCH, "intervals")

summary(unjoin(joined.ovenCH))

\dontrun{

## suppose the 5-year ovenbird covariates include a column for weight
## (here generated as random numbers)
for (i in 1:5) covariates(ovenCH[[i]])$wt <- runif(nrow(ovenCH[[i]]))
## construct single-session version of data for openCR
## identify 'wt' as varying across years
ovenCHj <- join(ovenCH, timevaryingcov = 'wt')
head(covariates(ovenCHj))
timevaryingcov(ovenCHj)
## Use example: openCR.fit(ovenCHj, model = p~wt)

}


}

\keyword{ manip }
