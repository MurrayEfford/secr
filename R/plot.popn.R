## plot.popn.R

## 2018-02-05

popIDsplit <- function (pop) {
    if (ms(pop)) {
        nsess <- length(pop)
        ID <-  unique(as.character(unlist(sapply(pop, rownames))))
        out <- array(dim = c(length(ID), nsess, 2), dimnames = list(ID, names(pop), c('x','y')))
        for (i in 1:nsess) {
            out[rownames(pop[[i]]),i,] <- unlist(pop[[i]])
        }
        out
    }
    else {
        stop("requires multisession pop object")    
    }
}

plot.popn <- function (x, add = FALSE, frame = TRUE, circles = NULL, collapse = FALSE, seqcol = NULL, ...) {
    if (ms(x)) {
        nsess <- length(x)
        temp <- do.call(rbind, lapply(x, function(y) attr(y,'boundingbox')))
        vertices <- apply(temp,2,range)
        if (collapse) {
            if (!add)  {
                MASS::eqscplot (0,0, xlab='', ylab='', xlim=vertices[,1],
                          ylim = vertices[,2], type='n', axes = FALSE)
            }
            if (frame) {    ## 2019-05-31, 2022-02-12
                if (!is.null(attr(x[[1]],'polygon'))) {
                    poly <- attr(x[[1]],'polygon')
                    if (inherits(poly, "SpatialPolygons")) {
                        poly <- st_as_sfc(poly)
                    }
                    if (inherits(poly, "sfc"))
                        plot(poly, add = TRUE)
                    else
                        polygon (poly)
                }
                else
                    polygon (attr(x[[1]], 'boundingbox'))   
            }
            sid <- popIDsplit(x)
            apply(sid, 1 ,lines, ...)
            if (!is.null(seqcol)) {
                nfc <- length(seqcol)
                if (nfc==1) seqcol <- rep(seqcol, nsess)
                else if (nfc==2) seqcol <- c(seqcol[1], rep(seqcol[2], nsess-1))
                else if (nfc != nsess) stop ("incorrect seqcol")
                first <- apply(!is.na(sid[,,1]), 1, match, x=TRUE)
                for (i in 1:nrow(sid)) {
                    x <- sid[i,first[i]:nsess,1]
                    y <- sid[i,first[i]:nsess,2]
                    points(x,y, pch=21, bg = seqcol[1:length(x)])
                }
            }
        }
        else {
            ## force shared frame
            for (i in 1:length(x)) attr(x,'boundingbox') <- vertices
            lapply (x, plot, add, frame, circles, ...)
        }
        invisible()
    }
    else {
        vertices <- attr(x,'boundingbox')
        
        if (!add)
        {
            if (frame)
                MASS::eqscplot (x$x, x$y, xlab='', ylab='', xlim=range(vertices$x),
                          ylim=range(vertices$y), type='n', axes = FALSE, ...)
            else
                MASS::eqscplot (x$x, x$y, xlab='', ylab='', type='n', axes = FALSE,
                          ...)
        }
        if (is.null(circles) | (nrow(x) == 0))    ## second condition 2011-09-14
            points (x$x, x$y, ...)
        else {
            if (length(circles) == 1)
                circles <- rep(circles, nrow(x))
            symbols (x$x, x$y, circles = circles, inches = FALSE,
                     add = TRUE, ...)
        }
        if (frame) {
            if (!is.null(attr(x,'polygon'))) {
                poly <- attr(x,'polygon')
                # we are assuming that the polygon attribute may be either
                # a 2-column matrix or a SpatialPolygons/sfc object
                if (inherits(poly, "SpatialPolygons")) {
                    poly <- st_as_sfc(poly)
                }
                if (inherits(poly, 'sfc')) {
                    plot(poly, add = TRUE)
                }
                else {
                    polygon (poly)
                }
            }
            else
                polygon (vertices)
        }
    }
}

