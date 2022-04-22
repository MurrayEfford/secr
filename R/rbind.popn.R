#############################################################################
## package 'secr'
## rbind.popn.R
## 2017-07-25 moved from methods.R
## 2017-07-25 now S3method; does not accept list input
#############################################################################

rbind.popn <- function (..., renumber = TRUE) {
    ## combine 2 or more popn objects
    ## ... may NOT be a single list object from 2017-07-25
    
    allargs <- list(...)
    
    if (any(sapply(allargs, inherits, "list")))
        stop ("rbind.popn no longer accepts list input; use do.call")
    
    ## check input
    check <- function (x) {
        if (!is(x,'popn'))
            stop ("all arguments must be 'popn' objects")
        if (is.null(covariates(x)) != is.null(covariates(allargs[[1]]) ))
            stop ("covariates must be provided for all or none")
    }
    sapply (allargs, check)
    
    ## construct output; default is hierarchical make.row.names
    animals <- rbind.data.frame(..., make.row.names = TRUE)
    
    ## row names
    if (renumber) {
        row.names(animals) <- 1:nrow(animals)
    }
    else {
        ## use original if unique, otherwise default to hierarchical
        an <- unlist(sapply(allargs, row.names, simplify = FALSE))
        if (any(duplicated(an))) 
            warning("renumbering popn to avoid duplicate row names")
        else
            row.names(animals) <- an
    }
    
    ## attributes
    class(animals) <- c('popn', 'data.frame')
    attr(animals, 'Ndist') <- 'user'
    attr(animals, 'model2D') <- attr(allargs[[1]], 'model2D')
    xl <- range(sapply(allargs, function(x) attr(x,'boundingbox')$x))
    yl <- range(sapply(allargs, function(x) attr(x,'boundingbox')$y))
    attr(animals, 'boundingbox') <- expand.grid (x=xl,y=yl)[c(1,3,4,2),]
    if (!is.null(covariates(allargs[[1]]))) {
        cov <- lapply(allargs, function(x) covariates(x))
        covariates(animals) <- do.call(rbind, cov)
        row.names(covariates(animals)) <- rownames(animals)
    }
    
    animals
}
###############################################################################
