\name{plot.mask}
\alias{plot.mask}
\alias{plot.Dsurface}
\alias{plot.Rsurface}
\alias{spotHeight}

\title{ Plot Habitat Mask, Density or Resource Surface }

\description{
  Plot a habitat mask either as points or as an \code{image}
  plot. Colours maybe used to show the value of one mask covariate.
}

\usage{
\method{plot}{mask}(x, border = 20, add = FALSE, covariate = NULL, axes = FALSE,
    dots = TRUE, col = "grey", breaks = 10, meshcol = NA, ppoly = TRUE,
    polycol = "red", legend = TRUE, \dots)

\method{plot}{Dsurface}(x, covariate, group = NULL, plottype =
    "shaded", scale = 1, ...)

\method{plot}{Rsurface}(x, covariate = "Resource", plottype =
    "shaded", scale = 1, ...)

spotHeight (object, prefix = NULL, dec = 2, point = FALSE, text = TRUE,
    sep = ", ",  session = 1, scale = 1, ...)

}

\arguments{
  \item{x, object}{ mask or Dsurface object }
  \item{border}{ width of blank display border (metres) }
  \item{add}{ logical for adding mask points to an existing plot }
  \item{covariate}{ name (as character string in quotes) or column
        number of a covariate to use for colouring }
  \item{axes}{ logical for plotting axes }
  \item{dots}{ logical for plotting mask points as dots, rather than as square pixels }
  \item{col}{ colour(s) to use for plotting }
  \item{breaks}{ an integer or a numeric vector -- see \code{\link{cut}} }
  \item{meshcol}{ colour for pixel borders (NA for none) }
  \item{ppoly}{ logical for whether the bounding polygon should be
    plotted (if `poly' specified) }
  \item{polycol}{ colour for outline of polygon (\code{ppoly} = TRUE) }
  \item{legend}{ logical; if TRUE a legend is plotted }   
  \item{\dots}{ other arguments passed to \code{eqscplot} (in the case
  of \code{plot.mask}), \code{plot.mask} (in the case of
  \code{plot.Dsurface} and \code{plot.Rsurface}), and \code{points} or
  \code{text} (in the case of \code{spotHeight}) }
  \item{group}{ group for which plot required, if more than 1 }
  \item{plottype}{ character string c("dots", "shaded", "contour",
    "persp") }
  \item{scale}{numeric multiplier for density or other numeric covariate
     (see \code{\link{Dsurface}})}
  \item{prefix}{character vector for name(s) of covariate(s) to retrieve}
  \item{dec}{ number of decimal places for rounding density }
  \item{point}{logical for whether to plot point}
  \item{text}{logical for whether to place density label on plot}
  \item{sep}{character separator for elements if length(prefix)>1}
  \item{session}{session number or identifier}
}

\details{

The argument \code{dots} of \code{plot.mask} selects between two
distinct types of plot (dots and shaded (coloured) pixels).

\code{plot.Dsurface} and \code{plot.Rsurface} offer contour and
perspective plots in addition to the options in \code{plot.mask}. It may
take some experimentation to get what you want - see
\code{\link[graphics]{contour}} and \code{\link[graphics]{persp}}.

For plot.Dsurface the default value of `covariate' is `D' unless the Dsurface has a `parameter' attribute of `noneuc',

If using a covariate or Dsurface or Rsurface to colour dots or pixels, the
\code{col} argument should be a colour vector of length equal to the
number of levels (the default palette from 2.9.0 is \code{terrain.colors}, and this
palette will also be used whenever there are too few levels in the
palette provided; see Notes for more on palettes). Border lines around
pixels are drawn in `meshcol'. Set this to NA to eliminate pixel
borders.

If a \code{covariate} is specified in a call to \code{plot.Dsurface} then
that covariate will be plotted instead of density. This is a handy way
to contour a covariate (contouring is not available in \code{plot.mask}).

If `breaks' is an integer then the range of the covariate is divided
into this number of equal intervals. Alternatively, `breaks' may be a
vector of break points (length one more than the number of
intervals). This gives more control and often `prettier' 

\code{spotHeight} may be used to interrogate a plot produced with
\code{plot.Dsurface} or \code{plot.Rsurface}, or by \code{plot.mask} if
the mask has covariates. \code{prefix} defaults to `density.' for
Dsurface objects and to `' (all covariates) for mask objects. The
predicted density or covariate at the nearest point is returned when the
user clicks on the plot. Multiple values may be displayed (e.g.,
\code{prefix = c("lcl","ucl")} if Dsurface includes confidence
limits). Click outside the mask or hit the Esc key to
end. \code{spotHeight} deals with one session at a time.

Legend plotting is enabled only when a covariate is specified. It uses
\code{\link{legend}} when \code{dots = TRUE} and
\code{\link{strip.legend}} otherwise.

}

\value{

  If \code{covariate} is specified and \code{plottype = "shaded"} then
  \code{plot.mask} invisibly returns a character vector of the intervals
  defined by `breaks' (useful for plotting a legend).

  If \code{plottype = "persp"} then \code{plot.mask} invisibly returns a the
  perspective matrix that may be used to add to the plot with
  \code{\link[grDevices]{trans3d}}.
  
  \code{spotHeight} invisibly returns a dataframe of the extracted
  values and their coordinates.
  
}

\seealso{
  \code{\link{colours}},
  \code{\link{mask}},
  \code{\link{Dsurface}},
  \code{\link{rectangularMask}},
  \code{\link[graphics]{contour}}
  \code{\link[graphics]{persp}}
  \code{\link{strip.legend}}
}

\note{

  \code{plot.mask()} acquired the argument `legend' in version 2.9.0,
  and other changes (e.g., \code{breaks = 10}) may alter the output.
  
  Contouring requires a rectangular grid; if a Dsurface is not
rectangular then plot.Dsurface with \code{plottype = "contour"} triggers a call to
\code{\link{rectangularMask}}.

The colour palettes \code{topo.colors}, \code{heat.colors} and
\code{terrain.colors} may be viewed with the \code{demo.pal} function in
the Examples code of their help page \link[grDevices]{palettes}.

The package \pkg{RColorBrewer} is a good source of palettes. Try
\code{display.brewer.all()} and e.g., \code{col = brewer.pal(7, "YlGn")}.

}

\examples{

# simple

temptrap <- make.grid()
tempmask <- make.mask(temptrap)
plot (tempmask)

\dontrun{

## restrict to points over an arbitrary detection threshold,
## add covariate, plot image and overlay traps

tempmask <- subset(tempmask, pdot(tempmask, temptrap,
    noccasions = 5)>0.001)
covariates (tempmask) <- data.frame(circle = 
    exp(-(tempmask$x^2 + tempmask$y^2)/10000) )
plot (tempmask, covariate = "circle", dots = FALSE, axes = TRUE, 
    add = TRUE, breaks = 8, col = terrain.colors(8), mesh = NA)
plot (temptrap, add = TRUE)

## add a legend
par(cex = 0.9)
covrange <- range(covariates(tempmask)$circle)
step <- diff(covrange)/8
colourlev <- terrain.colors(9)
zlev <- format(round(seq(covrange[1],covrange[2],step),2))
legend (x = "topright", fill = colourlev, legend = zlev, 
    y.intersp = 0.8, title = "Covariate")

title("Colour mask points with p.(X) > 0.001")
mtext(side=3,line=-1, "g0 = 0.2, sigma = 20, nocc = 5")

## Waitarere possum density surface extrapolated across region

regionmask <- make.mask(traps(possumCH), buffer = 1000, spacing = 10,
    poly = possumremovalarea)
dts <- distancetotrap(regionmask, possumarea)
covariates(regionmask) <- data.frame(d.to.shore = dts)
shorePossums <- predictDsurface(possum.model.Ds, regionmask)

## plot as coloured pixels with white lines
colourlev <- terrain.colors(7)
plot(shorePossums, breaks = seq(0,3.5,0.5), plottype = "shaded",
    poly = FALSE, col = colourlev, mesh = NA)
plot(traps(possumCH), add = TRUE, detpar = list(col = "black"))
polygon(possumremovalarea)

## check some point densities
spotHeight(shorePossums, dec = 1, col = "black")

## add a legend
zlev <- format(seq(0,3,0.5), digits = 1)
legend (x = "topright", fill = colourlev, legend =
    paste(zlev,"--"), y.intersp = 1, title = "Density / ha")
}

}

\keyword{ hplot }

