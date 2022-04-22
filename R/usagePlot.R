usagePlot <- function(object, add = FALSE, occasions = NULL, col =
                       'black', fill = FALSE, scale = 2, metres = TRUE, rad = 5, ...) {

    if (ms(object))  {
        lapply(object, usagePlot, add = add, occasions = occasions, col = col,
               fill = fill, scale = scale, metres = metres, rad = rad, ...)
    }
    else {
        if (is.null(usage(object)))
            stop ("object does not have usage attribute")
        if (!add)
            plot(object, ...)

        nocc <- ncol(usage(object))
        if (is.null(occasions)) {
            dx  <- rep((cos((1:nocc) * 2 * pi / nocc) * rad), each=nrow(object))
            dy  <- rep((sin((1:nocc) * 2 * pi / nocc) * rad), each=nrow(object))
            xy <- cbind(rep(object$x, nocc) + dx, rep(object$y, nocc)-dy)
            radius <- as.numeric(sqrt(usage(object)) * scale)
        }
        else {
            xy <- object
            if (tolower(occasions[1]) == "all")
                occasions <- 1:nocc
            usge <- apply(usage(object)[,occasions, drop = FALSE], 1, sum)
            radius <- sqrt(usge) * scale
        }

        # metres: use symbols with inches = FALSE
        if (metres) {
            fg <- col
            if (fill)
                bg <- fg
            else
                bg <- NULL
            radius[radius==0] <- NA  ## 2016-03-17
            symbols(xy, circles = radius, inches = FALSE, fg = fg, bg = bg, add = T)
        }
        else {
            pch <- ifelse (fill, 16, 1)
            points(xy, cex = radius, pch = pch, col = col)
        }
    }
}

bubble.legend <- function (legend, fill, col, 
                           height = 0.5, width = 0.06, 
                           text.cex = 0.9, xpd = TRUE, scale = 1, title = "", box = NA, box.col = par()$bg,
                           px = 0.95, py = 0.95, text.px = 0.06) {
    # cf strip.legend
    # check a plot exists?
    nleg <- length(legend)
    usr <- par()$usr
    if (length(py)!=2)
        py[2] <- py-0.2
        legendx <- (usr[2] - usr[1]) * px + usr[1]
        legendy <- (usr[4] - usr[3]) * seq(py[1],py[2], along.with=legend) + usr[3]
    
    legendx <- rep(legendx, nleg)
    symbols(legendx, legendy, circles = rev(legend)^0.5 * scale, inches = FALSE,
            add = TRUE, fg = col, bg = fill, xpd = TRUE)
    text(legendx + (usr[2] - usr[1])*text.px, legendy, rev(legend), adj = 0.5, xpd=TRUE)
}

####################################################################################
sightingPlot <- function (object, type = c("Detections", "Tu", "Tm", "Tn"),
                          add = FALSE, occasions = "ALL", mean = TRUE,
                          col = "black", fill = FALSE, scale = 2, metres = TRUE,
                          dropunused = TRUE, title = type,
                          legend = c(1,2,4,8), px = 0.95, py = 0.95, ...) {
    if (ms(object))  {
        lapply(object, sightingPlot, type = type, add = add, occasions = occasions,
               mean = mean, col = col, fill = fill, scale = scale, metres = metres, 
               dropunused = dropunused, title = title, legend = legend, px = px, 
               py = py, ...)
    }
    else {
        type <- match.arg(type)
        if (toupper(occasions)[1] != "ALL")
            object <- subset(object, occasions = occasions)
        if (dropunused & !is.null(usage(traps(object))))
            object <- suppressWarnings(subset(object, traps = 
                apply(usage(traps(object)),1,sum)>0))
        if (!add)
            plot(traps(object), ...)
        
        if (type == "Detections")
            tmp <- t(apply(object, 2:3, sum))
        else
            tmp <- attr(object, type)
        if (!is.null(usage(traps(object))))
            tmp[usage(traps(object))==0] <- NA
        
        ## aggregated in single bubble for each site
        if (mean)
            f <- apply(tmp, 1, mean, na.rm = TRUE)
        else
            f <- apply(tmp, 1, sum, na.rm = TRUE)
        radius <- sqrt(f) * scale
        
        xy <- traps(object)
        
        # metres: use symbols with inches = FALSE
        if (metres) {
            fg <- col
            if (fill)
                bg <- fg
            else
                bg <- NULL
            radius[radius==0] <- NA  ## 2016-03-17
            symbols(xy, circles = radius, inches = FALSE, fg = fg, bg = bg, add = T)
        }
        else {
            pch <- ifelse (fill, 16, 1)
            points(xy, cex = radius, pch = pch, col = col)
        }
        mtext(side = 3, title)
        
        ## add legend with various
        if (!is.null(legend)) {
            if (fill) fill <- col
            bubble.legend(legend=legend, fill=fill, col = col, scale = scale, px = px, py = py)
        }
        
        tr <- traps(object)
        covariates(tr) <- data.frame(f=f)
        invisible(tr)
    }
}
####################################################################################

