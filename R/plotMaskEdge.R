############################################################################################
## package 'secr'
## plotMaskEdge.R
############################################################################################

plotMaskEdge <- function (mask, plt = TRUE, add = FALSE, ...) {
    ## plots outer border of a mask, however irregular
    ## 2019-07-02 bug partly fixed by changing round(x) etc. to round(x,3) etc.
    ## 2020-08-26 now round(x,4) etc.
    draw <- function (x,y, vertical, OK) {
        if (!OK)
            rep(NA,4)
        else if (vertical)
            c(x, y+spc/2, x, y-spc/2)
        else
            c(x+spc/2, y, x-spc/2, y)
    }
    inmask <- function (x,y) {
        !is.na(match (paste(round(x,4),round(y,4),sep='.'), maskvector))
    }
    onecell <- function (xy) {
        x <- xy[1]
        y <- xy[2]
        cbind(
            draw (x-spc/2, y, TRUE, !inmask(x-spc, y)),
            draw (x+spc/2, y, TRUE, !inmask(x+spc, y)),
            draw (x, y-spc/2, FALSE, !inmask(x, y-spc)),
            draw (x, y+spc/2, FALSE, !inmask(x, y+spc))
        )
    }
    if (ms(mask)) {
        out <- lapply(mask, plotMaskEdge, plt=plt, add = add, ...)
        invisible(out)
    } else {
        if (plt & !add)
            plot(mask, dots = FALSE)
        spc <- spacing(mask)
        ## assumes mask has integer x,y resolution
        ## 2019-07-02 0.001 resolution
        ## 2020-08-26 0.0001 resolution
        maskvector <- apply(round(mask,4),1,paste, collapse='.')
        coord <- apply(mask, 1, onecell)
        xy <- matrix(coord, nrow = 4)
        xy <- xy[,!apply(xy, 2, function(z) any(is.na(z)))] # drop NA
        if (plt) segments(xy[1,], xy[2,], xy[3,],xy[4,], ...)
        invisible(xy)
    }
}

#########################################
