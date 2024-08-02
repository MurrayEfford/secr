# MCgof methods

plot.MCgof <- function(x, overlay = NULL, ...) {
    main <- c('individual-detector','individual', 'detector')
    onestat <- function (xy, ynum) {
        lim <- range(xy[1:2,])
        lim <- lim + c(-1, +1) * diff(lim)/7
        MASS::eqscplot(xy['Tsim',], xy['Tobs',], 
                       xlab = 'simulated', ylab = 'observed',
                       xlim = lim, ylim = lim)
        
        abline(0,1, col = 'red')
        mtext(side=3, paste(main[ynum], " p =", round(mean(xy[3,]),3)))
        
        # optional overlay of points from another MCgof or from scrgof
        if (!is.null(overlay)) {
            if (inherits(overlay, 'MCgof')) {
                xy <- overlay$all[[ynum]]
                points(xy['Tsim',], xy['Tobs',], ...)
            }
            else if (names(overlay)[1] == 'scrgof_pval') {
                points(overlay[[ynum+1]], ...)
            }
        }        
    }
    par(pty = 's', mfrow = c(2,2))
    mapply(onestat, x$all, 1:3)
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
