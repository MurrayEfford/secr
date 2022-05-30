# 2022-05-15

# inputs (D,sigma) vector or fitted model or dataframe from predict or list of preceding
# return k = sigma / 100 * sqrt(D) assuming sigma in metres, D in animals/hectare

kfn <- function (object) { 
    if (is.matrix(object)) {
        out <- t(apply(object,1,kfn))
        dimnames(out) <- list(rownames(object), c('D','sigma','k'))
        out
    }
    else if (is.numeric(object)) {
        D <- object[1]
        sigma <- object[2]
        c(D = D, sigma = sigma, k = sigma/100 * sqrt(D))    
    }
    else if (is.data.frame(object)) {
        if (!all(c('sigma','D') %in% rownames(object))) stop()
        sigma <- object['sigma', 'estimate']
        D <- object['D','estimate']
        c(D = D, sigma = sigma, k = sigma/100 * sqrt(D))    
    }
    else {
        if (inherits(object, 'secr')) {
            if (object$detectfn %in% c('HN','HHN')) warning("fitted with non-normal detectfn")
            object <- predict(object)
            kfn(object)
        }
        else {
            if (is.data.frame(object[[1]])) {
                t(sapply(object, kfn))
            }
            else {
                # assume list of lists of data.frames
                lapply(object, function(x) t(sapply(x, kfn)))
            }
        }
    }
}

# kfn(secrdemo.0)
# [1] 0.6874
