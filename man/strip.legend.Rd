\name{strip.legend}
\alias{strip.legend}

\title{Colour Strip Legend}

\description{

  This function is used with shaded plots to display a legend.

}

\usage{

strip.legend(xy, legend, col, legendtype = c("breaks", "intervals", "other"),
  tileborder = NA, height = 0.5, width = 0.06, inset = 0.06, text.offset = 0.02,
  text.cex = 0.9, xpd = TRUE, scale = 1, title = "", box = NA, box.col = par()$bg)

}

\arguments{
  \item{xy}{location of legend (see Details)}

\item{legend}{character vector (see Details)}

\item{col}{vector of colour values}

\item{legendtype}{character}

\item{tileborder}{colour of lines around each tile in the colour strip. Use
  NA for none.}

\item{height}{height of colour strip as a fraction of the plot dimensions}

\item{width}{width of colour strip as a fraction of the plot dimensions}

\item{inset}{spacing between legend and outside plot boundary, as a
fraction of the plot dimensions}

\item{text.offset}{spacing between colour strip and text, as a
fraction of the plot dimensions}

\item{text.cex}{size of text font}

\item{xpd}{logical, if TRUE the legend will use the margins of the plot}

\item{scale}{numeric; each value x will be displayed as scale * x}

\item{title}{text displayed above legend}

\item{box}{colour of frame, if framed, otherwise NA}

\item{box.col}{colour of background, if framed, otherwise ignored}

}

\details{
  
The location of the legend is determined by \code{xy} which may be one
of the character values "topright", "topleft", "bottomright",
"bottomleft", "right", "left", or the x-y coordinates (in
user units) of the top-left corner of the colour strip. Coordinates may
be given as a vector or a list, and the output from
\code{\link{locator}}(1) is suitable.

For more on colours, see notes in \code{\link{plot.mask}} and
\code{\link{colors}} and \code{\link{terrain.colors}}

If \code{legendtype = 'breaks'} then labels are placed at the class
boundaries; otherwise, the labels are centred vertically. If
\code{legendtype = 'breaks'} or  \code{legendtype = 'intervals'} then
numeric values are extracted from the input, otherwise the text strings
in \code{legend} are used as provided. 

The legend itself may be provided as a vector of values or as the class
labels output from \code{\link{plot.mask}}. Class labels are generated
by \code{\link{cut}} in the form `(0,20]', `(20,40]', etc. These are
parsed to construct either breaks (0,20,40,...) or intervals (`0-20',
`20-40',...) as requested in the \code{legendtype} argument.

\code{box} may also be TRUE/FALSE; if TRUE the foreground colour is used \code{par()$fg}.

}

\value{
  
  Invisibly returns a vector of user coordinates for the left, right,
  bottom and top of the colour strip.

}

\note{
  
From \pkg{secr} 2.9.0, the default behaviour of \code{plot.mask} is to
call \code{strip.legend} to display a legend in the top right of the
plot, labeled at breaks.

}

\seealso{

  \code{\link{plot.mask}}
  
}

\examples{


temptrap <- make.grid()
tempmask <- make.mask(temptrap)
covariates (tempmask) <- data.frame(circle = 
    exp(-(tempmask$x^2 + tempmask$y^2)/10000) )
tmpleg <- plot (tempmask, covariate = "circle", dots = FALSE, 
    breaks = 10, legend = FALSE)
strip.legend (xy = 'topright', col = terrain.colors(10),
    legend = tmpleg, title = "Test plot")

if (interactive()) {
    ## a custom axis using the returned values
    par(mar = c(2,2,2,6))
    plot (tempmask, covariate = "circle", dots = FALSE, 
        breaks = 10, legend = FALSE)
    b <- strip.legend (locator(1), col = terrain.colors(10),
        legendtype = "other", legend = "    ", title = "Test plot",
         height = 0.3, box = NA)
    axis(side = 4, pos = b[2]+5, at = seq(b[4], b[3], length = 3),
        lab = seq(0,1,0.5), las = 1, tck = -0.02)
    par(mar = c(5,4,4,2) + 0.1)  ## reset to default
}

}

