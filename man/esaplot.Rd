\name{esaPlot}
\alias{esaPlot}

\title{
Mask Buffer Diagnostic Plot
}

\description{

Plot effective sampling area (Borchers and Efford 2008) as a function of increasing buffer width.

\code{esaPlot} was previously called \code{esa.plot}.

}

\usage{

esaPlot (object, max.buffer = NULL, spacing = NULL, max.mask = NULL,
    detectfn, detectpar, noccasions, binomN = NULL, thin = 0.1,
    poly = NULL, poly.habitat = TRUE, session = 1, plt = TRUE, 
    type = c('density', 'esa', 'meanpdot', 'CVpdot'), n = 1, add = FALSE, 
    overlay = TRUE, conditional = FALSE, ...)

}

\arguments{
  \item{object}{\code{traps} object or \code{secr} object output from \code{secr.fit} }
  \item{max.buffer}{ maximum width of buffer in metres }
  \item{spacing}{ distance between mask points }
  \item{max.mask}{ \code{\link{mask}} object }
  \item{detectfn}{ integer code or character string for shape of
  detection function 0 = halfnormal etc. -- see \link{detectfn}}
  \item{detectpar}{ list of values for named parameters of detection function }
  \item{noccasions}{ number of sampling occasions}
  \item{binomN}{ integer code for discrete distribution (see
  \code{\link{secr.fit}}) }  
  \item{thin}{ proportion of mask points to retain in plot and output}
  \item{poly}{ matrix of two columns interpreted as the x and y
    coordinates of a bounding polygon (optional) }
  \item{poly.habitat}{logical as in \code{\link{make.mask}} }
  \item{session}{vector of session indices (used if \code{object} spans multiple sessions)}
  \item{plt}{ logical to plot results }
  \item{type}{ character, what to plot }
  \item{n}{ integer number of distinct individuals detected }
  \item{add}{ logical to add line to an existing plot }
  \item{overlay}{ logical; if TRUE then automatically \code{add = TRUE}
  for plots after the first }
  \item{conditional}{logical; if TRUE the reported mean and CV are conditional on detection 
      (see \code{\link{CVpdot}})}
  \item{\dots}{ graphical arguments passed to plot() and lines()}
}

\details{

  Effective sampling area (esa) is defined as the integral of net
  capture probability (\eqn{p.(\mathbf{X})}{p.(X)}) over a
  region. \code{esaPlot} shows the effect of increasing region size on
  the value of esa for fixed values of the detection parameters. The
  \code{max.buffer} or \code{max.mask} arguments establish the maximum
  extent of the region; points (cells) within this mask are sorted by
  their distance \eqn{d_k}{dk} from the nearest detector. esa(buffer) is
  defined as the cumulative sum of \eqn{cp.(\mathbf{X})}{c * p.(X)} for
  \eqn{d_k(\mathbf{X}) <= \mathrm{buffer}}{dk(X) <= buffer}, where \eqn{c}
  is the area associated with each cell.

  The default (\code{type = 'density'}) is to plot the reciprocal of esa
  multiplied by \code{n}; this is on a more familiar scale (the density
  scale) and hence is easier to interpret.
  
  Because \code{esaPlot} uses the criterion `distance to nearest
  detector', \code{max.mask} should be constructed to include all
  habitable cells within the desired maximum buffer and no others. This
  is achieved with \code{type = "trapbuffer"} in \code{make.mask}. It is
  a good idea to set the \code{spacing} argument of \code{make.mask}
  rather than relying on the default based on \code{nx}. Spacing may be
  small (e.g. sigma/10) and the buffer of \code{max.mask} may be quite
  large (e.g. 10 sigma), as computation is fast.
  
  Thinning serves to reduce redundancy in the plotted points, and (if
  the result is saved and printed) to generate more legible numerical
  output. Use \code{thin=1} to include all points.

  \code{esaPlot} calls the internal function \code{esaPlotsecr} when
  \code{object} is a fitted model. In this case \code{detectfn},
  \code{detectpar} and \code{noccasions} are inferred from
  \code{object}. 
}

\value{
  
  A dataframe with columns
  \item{buffer}{ buffer width }
  \item{esa}{ computed effective sampling area }
  \item{density}{ n/esa } 
  \item{pdot}{ \eqn{p.(\mathbf{X})}{p.(X)} }
  \item{pdotmin}{ cumulative minimum  (\eqn{p.(\mathbf{X})}{p.(X)}) }
  \item{meanpdot}{ expected pdot across mask (see \code{\link{CVpdot}})}
  \item{CVpdot}{ CV of pdot across mask (see \code{\link{CVpdot}})}
  
  If \code{plt = TRUE} the dataframe is returned invisibly.

}

\note{
The response of effective sampling area to buffer width is just one
possible mask diagnostic; it's fast, graphic, and often
sufficient. \code{\link{mask.check}} performs more intensive checks,
usually for a smaller number of buffer widths.

The old argument 'as.density' was superceded by 'type' in 3.1.7.

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

}

\seealso{ \code{\link{mask}}, \code{\link{pdot}}, \code{\link{CVpdot}},
  \code{\link{make.mask}}, \code{\link{mask.check}},
  \code{\link{Detection functions}} }

\examples{

\dontrun{

## with previously fitted model
esaPlot(secrdemo.0)

## from scratch
trps <- make.grid()
msk <- make.mask(trps, buffer = 200, spacing = 5, type = "trapbuffer")
detectpar <- list(g0 = 0.2, sigma = 25)
esaPlot(trps,,, msk, 0, detectpar, nocc = 10, col = "blue")
esaPlot(trps,,, msk, 0, detectpar, nocc = 5, col = "green",
    add = TRUE)

esaPlot(trps,,, msk, 0, detectpar, nocc = 5, thin = 0.002, plt = FALSE)

}

}

% Not exported
% 
% esaPlotsecr (object, max.buffer = NULL, max.mask = NULL,
%     thin = 0.1, poly = NULL, poly.habitat = TRUE, session = 1, plt = TRUE, 
%     type = "density", add = FALSE, overlay = TRUE, conditional = FALSE, ...)
% 
% }
% 
% \arguments{
% 
%   \item{object}{\code{secr} object output from \code{secr.fit} }
%   \item{max.buffer}{ maximum width of buffer in metres }
%   \item{max.mask}{ \code{\link{mask}} object }
%   \item{thin}{ proportion of mask points to retain in plot and output}
%   \item{poly.habitat}{logical as in \code{\link{make.mask}} }
%   \item{session}{vector of session indices (used if \code{object} spans multiple sessions)}
%   \item{plt}{ logical to plot results }
%   \item{type}{ character; see \code{\link{esaPlot}}}
%   \item{add}{ logical to add line to an existing plot }
%   \item{overlay}{ logical; if TRUE then automatically \code{add = TRUE}
%   for plots after the first }
%   \item{conditional}{logical; see \code{\link{esaPlot}}}
%   \item{\dots}{ other graphical arguments passed to plot() and lines() (esaPlotsecr)}
% 
% }
% 
% \details{
% 
%   \code{esaPlotsecr} provides a wrapper for \code{\link{esaPlot}}
%   that is called internally from esaPlot when it is presented with an
%   secr object. Arguments of \code{esaPlot} such as \code{detectfn} are
%   inferred from the fitted model.
% 
%   The optional `poly' argument of \code{esaPlotsecr} is a matrix of two columns 
%   interpreted as the x and y coordinates of a bounding polygon.
% 
%   If \code{max.mask} is not specified then a maximal mask of type
%   `trapbuffer' is constructed using \code{max.buffer} and the spacing of
%   the mask in \code{object}. In this case, if \code{max.buffer} is not
%   specified then it is set either to the width of the existing plot
%   (\code{add = TRUE}) or to 10 x sigma-hat from the fitted model in
%   \code{object} (\code{add = FALSE}).
% 
% }
