## secr
## 2021-07-12

extractMoves <- function (pop, plotn = NULL, add = FALSE, collapse = TRUE, maxradius = Inf, ...) {
    # If plotn is provided then the first plotn movements will be 
    # plotted for each successive pair of sessions; maxradius restricts 
    # initial locations to within that distance of the centre.
    npop <- length(pop)
    if (npop<2 || !inherits(pop, 'popn')) {
        stop ("requires multisession popn object")
    }
    centre <- apply(attributes(pop[[1]])$boundingbox, 2, mean)
    d <- vector('list', npop-1)
    for (i in 1:(npop-1)) {
        xy <- sweep(pop[[i]], MARGIN = 2, FUN = "-", STATS = centre)
        r1OK <- apply(xy^2,1,sum)^0.5 < maxradius
        pop0 <- subset(pop[[i]], r1OK)
        pop1 <- pop0[rownames(pop0) %in% rownames(pop[[i+1]]),]
        pop2 <- pop[[i+1]][rownames(pop[[i+1]]) %in% rownames(pop0),]
        if (!is.null(plotn)) {
            tempn <- min(nrow(pop1), plotn)
            if (!add && (i == 1 || !collapse)) plot(pop1[1:tempn,])
            arrows(
                pop1[1:tempn,]$x, 
                pop1[1:tempn,]$y, 
                pop2[1:tempn,]$x, 
                pop2[1:tempn,]$y,
                ...)
            points(pop1[1:tempn,], pch = 16, xpd = TRUE)
            points(pop2[1:tempn,], pch = 16, xpd = TRUE)
        }
        d[[i]] <- cbind(pop1, pop2, apply((pop1-pop2)^2,1,sum)^0.5)
        names(d[[i]]) <- c('x1','y1','x2','y2','d')
    }
    d
}
