# MCgof methods

plot.MCgof <- function(x, counts = 'all', overlay = NULL, ...) {
    onestat <- function (count) {
        xy <- x$all[[count]]    # dimensions c(Tsim, Tobs, p), 1:nreplicates)
        lim <- range(xy[1:2,])
        lim <- lim + c(-1, +1) * diff(lim)/7
        MASS::eqscplot(xy['Tsim',], xy['Tobs',], 
                       xlab = 'simulated', ylab = 'observed',
                       xlim = lim, ylim = lim)
        
        abline(0,1, col = 'red')
        mtext(side=3, paste(main[count], " p =", round(mean(xy[3,]),3)))
        
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
    main <- c(yik = 'individual-detector', yi = 'individual', yk = 'detector')
    if (tolower(counts) == 'all') counts <- names(x$all)
    junk <- lapply(counts, onestat)
    invisible()
}
################################################################################

summary.MCgof <- function(object, ...) {
    temp <- object[c('nsim', 'means','proctime')]
    class(temp) <- 'summary.MCgof'
    temp
}
################################################################################

print.summary.MCgof <- function (x, ...) {
    cat (paste0('nsim = ', x$nsim, ', proctime = ', round(x$proctime, 2), ' seconds\n'))
    print(x$means)
}
################################################################################
