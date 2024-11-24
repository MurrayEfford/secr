############################################################################################
## package 'secr'
## split.mask.R
## last changed 2018-12-18
## 2024-11-08 check for NA in f
############################################################################################

split.mask <- function (x, f, drop = FALSE, clusters = NULL, ...) {
    if (ms(x))
        stop ("split not suitable for multi-session mask")
    
    if (!is.null(clusters)) {
        if (!missing(f))
            stop ("provide only one of f or clusters")
        if (!ms(clusters) || ! inherits(clusters, 'traps'))
            stop ("clusters should be list of traps objects")
        buff <- max(distancetotrap(x, rbind(clusters)))
        nclusters <- length(clusters)
        onecluster <- function (cl) {
            subset(x, distancetotrap(x, clusters[[cl]]) <= buff) 
        }
        out <- lapply(1:nclusters, onecluster)
        names(out) <- names(clusters)
    } 
    else {
        f <- if (drop) factor(f) else as.factor(f)
        if (length(f) != nrow(x))
            stop ("length of f should match number of rows in mask")
        if (any(is.na(f))) {
            stop ("f should not contain NA values")
        }        
        out <- list()
        for (i in levels(f)) {
            out[[i]] <- subset (x, subset = (f == i), ...)
        }
        names(out) <- levels(f)
    }
    class(out) <- c('mask','list')
    out
}
