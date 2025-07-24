##############################################################################
## package 'secr'
## plot.traps.R
## 2017-01-26 moved from methods.R
## 2022-10-19 add frame argument
##############################################################################

plot.traps <- function(x,
    border = 100,
    label = FALSE,
    offset = c(6,6),
    add = FALSE,
    hidetr = FALSE,
    detpar=list(),
    txtpar=list(),
    bg = 'white',
    gridlines = !add,
    gridspace = 100,
    gridcol = 'grey',
    markused = FALSE,
    markvarying = FALSE,
    markvertices = FALSE,
    labelclusters = FALSE,
    frame = NULL,
    ... )
{
    #### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12
    if (ms(x)) {
        lapply(x, plot.traps,
            border, label, offset,
            add, hidetr, detpar, txtpar,
            bg, gridlines, gridspace, gridcol,
            markused, markvarying, markvertices, labelclusters, 
            frame, ...)
    }
    else {
        trappar <- list(...)
        buff <- c(-border,+border)
        if (add) {
            xl <-  par()$usr[1:2]  
            yl <-  par()$usr[3:4]  
        } else {
            xl <- range(x$x)+buff
            yl <- range(x$y)+buff
        }
        offsety <- ifelse (length(offset)==2, offset[2], offset[1])
        if (all(detector(x) %in% c('polygon','polygonX')))
            dcol <- 'white'
        else 
            dcol <- 'red'
        detpar <- secr_replacedefaults (list(col=dcol, pch=3, cex=0.8, fg='red', bg='white'), detpar)
        txtpar <- secr_replacedefaults (list(col='blue', cex=0.7), txtpar)
        if (is.logical(markvertices))
            markvertices <- markvertices * 2  ## 0 or 2
        
        if (!is.null(usage(x))) {
            used <- apply(attr(x,'usage'),1,function(z) any(z>0))
            varying <- used * apply(attr(x,'usage'),1,function(z) any(z==0))
        }
        else {
            used  <- rep(TRUE, nrow(x))
            varying <- rep(FALSE, nrow(x))
        }
        initialpar <- par(detpar)
        on.exit(par(initialpar))
        
        if (!add) {
            par(bg=bg)
            ## axes = FALSE blocks bty = 'o' 2011-05-08
            if (!is.null(frame)) {
                plot(boundarytoSF(frame), border = 'black')
            }
            else {
                MASS::eqscplot (x$x, x$y, xlim=range(x$x)+buff, ylim=range(x$y)+buff,
                    xlab='', ylab='', type='n', axes=F, ...)
            }
            if (!is.null(trappar$bty)) {
                if (trappar$bty=='o') rect(xl[1],yl[1],xl[2],yl[2], border = 'black')
            }
        }
        plotvertices <- function (df) {
            if (markvertices == 1)
                i <- c(1, nrow(df))
            else
                i <- 1:nrow(df)
            points(df$x[i], df$y[i], pch = detpar$pch, bg = bg, col=detpar$col)
        }
        if (gridlines) {
            strtx <- gridspace * floor(xl[1]/gridspace)
            strty <- gridspace * floor(yl[1]/gridspace)
            finx  <- gridspace * (floor(xl[2]/gridspace) + 1)
            finy  <- gridspace * (floor(yl[2]/gridspace) + 1)
            for (xi in seq(strtx, finx, gridspace)) {
                if (xi>=xl[1] & xi<= xl[2]) {
                    x1 <- x2 <- xi
                    y1 <- max(strty, yl[1])
                    y2 <- min(finy, yl[2])
                    segments(x1, y1, x2, y2, col = gridcol)
                }
            }
            for (yi in seq(strty, finy, gridspace)) {
                if (yi>=yl[1] & yi<=yl[2]) {
                    y1 <- y2 <- yi
                    x1 <- max(strtx, xl[1])
                    x2 <- min(finx, xl[2])
                    segments(x1, y1, x2, y2, col = gridcol)
                }
            }
        }
        
        if (!hidetr) {
            if (all(detector(x) %in% c('polygon','polygonX'))) {
                templist <- split (x, levels(polyID(x)), prefix='')
                lapply(templist, function (y) {
                    polygon (y$x, y$y, col = detpar$col, 
                        density = if (is.na(detpar$col)) 0 else NA, 
                        border = detpar$fg)
                    })
                if (markvertices > 0) {
                    lapply(templist, plotvertices)
                }
                if (label) for (k in 1:length(templist)) {
                    if (all(detector(x) %in% c('polygon','polygonX'))) {
                        msk <- suppressWarnings(make.mask(templist[[k]], buffer = 0, poly = templist[[k]], nx = 32))
                        xbar <- mean(msk$x)
                        ybar <- mean(msk$y)
                    }
                    else {
                        xbar <- mean(range(templist[[k]]$x))
                        ybar <- mean(range(templist[[k]]$y))
                    }
                    text (xbar+offset[1], ybar+offsety, names(templist)[k])
                }
            }
            else if (all(detector(x) %in% c('transect','transectX'))) {
                templist <- split (x, levels(transectID(x)), prefix='')
                lapply(templist, function (df) lines (df$x, df$y, col=detpar$col))
                if (markvertices > 0) {
                    lapply(templist, plotvertices)
                }
                if (label) for (k in 1:length(templist)) {
                    xbar <- mean(range(templist[[k]]$x))
                    ybar <- mean(range(templist[[k]]$y))
                    text (xbar+offset[1], ybar+offsety, names(templist)[k])
                }
            }
            else {
                points (x$x, x$y)
                if (markused) {
                    points (x$x[used], x$y[used], pch = 1, cex = 0.8)
                }
                if (markvarying & any(varying)) {
                    points (x$x[varying], x$y[varying], pch = 16, cex = 0.8)
                }
            }
            par(txtpar)
            if (label & !(all(detector(x) %in% .localstuff$polydetectors))) {
                text (x$x+offset[1], x$y+offsety, rownames(x))
            }
            if (labelclusters & !all(detector(x) %in% .localstuff$polydetectors)) {
                if (is.null(clusterID(x)) | is.null(clustertrap(x)))
                    stop ("require clustered traps to label with clusterID")
                cl1 <- clustertrap(x) == 1
                text (x$x[cl1]+offset[1], x$y[cl1]+offsety, clusterID(x)[cl1])
            }
            #par(initialpar)   # restore
        }
        invisible()
    }
}
###############################################################################
