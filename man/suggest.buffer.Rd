\name{suggest.buffer}
\alias{suggest.buffer}
\alias{bias.D}
\alias{bufferBiasCheck}

\title{
Mask Buffer Width
}

\description{

Determines a suitable buffer width for an integration \link{mask}. The
`buffer' in question defines a concave polygon around a detector array
constructed using \code{make.mask} with \code{type = "trapbuffer"}.  The
method relies on an approximation to the bias of maximum likelihood
density estimates (M. Efford unpubl).

\code{bufferBiasCheck} is the post hoc test implemented in 
\code{\link{secr.fit}} when 'biasLimit' is not NA.

}

\usage{
suggest.buffer(object, detectfn = NULL, detectpar = NULL, 
    noccasions = NULL, ignoreusage = FALSE, ncores = NULL, RBtarget = 0.001,
    interval = NULL, binomN = NULL, ...)

bias.D (buffer, traps, detectfn, detectpar, noccasions, binomN = NULL,
    control = NULL)

bufferBiasCheck(object, buffer, biasLimit)

}

\arguments{
  \item{object}{ single-session `secr', `traps' or `capthist' object }
  \item{detectfn}{ integer code or character string for shape of
  detection function 0 = halfnormal etc. -- see \link{detectfn}}
  \item{detectpar}{ list of values for named parameters of detection
  function -- see \link{detectpar} }
  \item{noccasions}{ number of sampling occasions}
  \item{ignoreusage}{logical for whether to discard usage information from
    \code{traps(capthist)}}
  \item{ncores}{ integer number of threads to use for parallel processing}
  \item{RBtarget}{ numeric target for relative bias of density estimate }
  \item{interval}{ a vector containing the end-points of the interval to
    be searched }
  \item{binomN}{ integer code for distribution of counts (see
  \code{\link{secr.fit}}) }
  \item{...}{ other argument(s) passed to \code{bias.D} }
  \item{buffer}{ vector of buffer widths }
  \item{traps}{ `traps' object }
  \item{control}{ list of mostly obscure numerical settings (see Details) }
  \item{biasLimit}{ numeric threshold for predicted relative bias due to
  buffer being too small }
}
\details{

The basic input style of \code{suggest.buffer} uses a `traps' object and
a detection model specified by `detectpar', `detectfn' and `noccasions',
plus a target relative bias (RB). A numerical search is conducted for
the buffer width that is predicted to deliver the requested RB. If
\code{interval} is omitted it defaults to (1, 100S) where S is the
spatial scale of the detection function (usually
\code{detectpar$sigma}). An error is reported if the required buffer
width is not within \code{interval}. This often happens with
heavy-tailed detection functions (e.g., hazard-rate): choose another
function, a larger \code{RBtarget} or a wider \code{interval}.

Setting \code{ncores = NULL} uses the existing value from the environment variable 
RCPP_PARALLEL_NUM_THREADS (see \code{\link{setNumThreads}}).

Convenient alternative input styles are --
\itemize{
  
\item \code{secr} object containing a fitted model. Values of `traps',
`detectpar', `detectfn' and `noccasions' are extracted from
\code{object} and any values supplied for these arguments are ignored.

\item \code{capthist} object containing raw data. If \code{detectpar} is
not supplied then \code{\link{autoini}} is used to get `quick and dirty'
values of \code{g0} and \code{sigma} for a halfnormal detection
function. \code{noccasions} is ignored. \code{autoini} tends to
underestimate \code{sigma}, and the resulting buffer also tends to be
too small.

}

\code{bias.D} is called internally by \code{suggest.buffer}.
% gpclib suppressed 2012-11-03
% The package \pkg{gpclib} may be used for more accurate estimates of the length of buffer
% contours (this does not appear to be critical). Some components of
% \code{control} are specific to this part of the algorithm
% (\code{ntheta}, \code{ninterp}, \code{maxinterp}).
}
\value{

  \code{suggest.buffer} returns a scalar value for the suggested buffer
  width in metres, or a vector of such values in the case of a
  multi-session \code{object}.

  \code{bias.D} returns a dataframe with columns \code{buffer} and \code{RB.D}
  (approximate bias of density estimate using finite buffer width,
  relative to estimate with infinite buffer).
  
  \code{bufferBiasCheck} returns a logical value, FALSE if the predicted bias 
  exceeds biasLimit.
}

\note{
  
The algorithm in \code{bias.D} uses one-dimensional numerical
integration of a polar approximation to site-specific detection
probability. This uses a further 3-part linear approximation for the
length of contours of distance-to-nearest-detector (\eqn{r}) as a
function of \eqn{r}.

The approximation seems to work well for a compact detector array, but
it should not be taken as an estimate of the bias for any other purpose:
do \emph{not} report \code{RB.D} as "the relative bias of the density
estimate". \code{RB.D} addresses only the effect of using a finite
buffer. The effect of buffer width on final estimates should be checked
with \code{\link{mask.check}}.

The default buffer type in \code{make.mask}, and hence in
\code{secr.fit}, is `traprect', not `trapbuffer', but a buffer width
that is adequate for `trapbuffer' is always adequate for `traprect'.

\code{control} contains various settings of little interest to the
user.

The potential components of \code{control} are --

\describe{
  \item{\code{method = 1 }}{ code for method of modelling \emph{p.(X)} as a
    function of buffer (\emph{q(r)})}
  
  \item{\code{bfactor = 20 }}{ \emph{q(r)} vs \emph{p.(X)} calibration
    mask buffer width in multiples of trap spacing }

  \item{\code{masksample = 1000 }}{ maximum number of points sampled from
    calibration mask }
  
  \item{\code{spline.df = 10 }}{ effective degrees of freedom for
    \code{smooth.spline} }
    
    \item{\code{ncores = NULL}}{ integer number of cores}
  
%  \item{\code{use.gpclib = FALSE }}{ logical to use \pkg{gpclib} if
%    installed } 
%  \item{\code{ntheta = 60 } }{ integer value for smoothness of contours
%    (straight-line approximation to arc lengths)}  
%  \item{\code{ninterp = 5 } }{ number of points to interpolate between
%  trapspacing/2 and trapspacing/2 ^0.5 on the contour-length vs buffer
%  curve }
}

}

\seealso{

\code{\link{mask}}, \code{\link{make.mask}}, \code{\link{mask.check}}, \code{\link{esaPlot}}

}

\examples{

\dontrun{

temptraps <- make.grid()
detpar <-  list(g0 = 0.2, sigma = 25)
suggest.buffer(temptraps, "halfnormal", detpar, 5)

suggest.buffer(secrdemo.0)

suggest.buffer(ovenCH[[1]])

RB <- bias.D(50:150, temptraps, "halfnormal", detpar, 5)
plot(RB)

detpar <-  list(g0 = 0.2, sigma = 25, z=5)
RB <- bias.D(50:150, temptraps, "hazard rate", detpar, 5)
lines(RB)

## compare to esa plot
esaPlot (temptraps, max.buffer = 150, spacing = 4, detectfn = 0,
  detectpar = detpar, noccasions = 5, type = "density")

## compare detection histories and fitted model as input
suggest.buffer(captdata)
suggest.buffer(secrdemo.0)

}

}

