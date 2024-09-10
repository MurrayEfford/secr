# MCgof methods

plot.MCgof <- function(x, counts = 'all', overlay = NULL, maxT = NULL, main = NULL, cex = 0.9, ...) {
    onestat <- function (count) {
        xy <- x$all[[count]]    # dimensions c(Tsim, Tobs, p), 1:nreplicates)
        if (!is.null(maxT)) xy[xy>maxT] <- NA
        lim <- range(xy[1:2,], na.rm = TRUE)
        lim <- lim + c(-1, +1) * diff(lim)/7
        MASS::eqscplot(xy['Tsim',], xy['Tobs',], 
                       xlab = 'simulated', ylab = 'observed',
                       xlim = lim, ylim = lim)
        
        abline(0,1, col = 'red')
        if (any(main != "")) {
            mtext(side=3, line=0.4, paste(main[count], " p =", round(summ['p',count],3)), cex = cex)
        }
        
        # optional overlay of points from another MCgof or from scrgof
        if (!is.null(overlay)) {
            if (inherits(overlay, 'MCgof')) {
                xy <- overlay$all[[count]]
                points(xy['Tsim',], xy['Tobs',], ...)
            }
            else if (names(overlay)[1] == 'scrgof_pval') {
                scrgofnames <- c(yik = 'gof_ik', yi = 'gof_i', yk = 'gof_j')
                points(overlay[[scrgofnames[count]]], ...)
            }
        }        
    }
    if (is.null(main)) {
        # default vector of potential labels
        main <- c(yik = 'yik individual x detector', 
                  yi = 'yi individual', 
                  yk = 'yk detector')
    }
    opar <- par(cex = cex)
    on.exit(par(opar))
    if (tolower(counts[1]) == 'all') counts <- names(x$all)
    summ <- summary(x)
    junk <- lapply(counts, onestat)
    invisible()
}
################################################################################

summary.MCgof <- function(object, ...) {
    summ <- sapply(object$all, apply, 1, median, na.rm = TRUE)
    summ['p',] <- sapply(object$all, function(x) mean(x['p',], na.rm = TRUE))
    summ <- rbind(summ, validn = sapply(object$all, function(x) sum(!is.na(x['p',]))))
    class(summ) <- 'summary.MCgof'
    summ
}
################################################################################

print.summary.MCgof <- function (x, ...) {
    class(x) <- NULL
    print(x)
}
################################################################################

print.MCgof <- function (x, ...) {
    print(summary(x))
}
################################################################################
