#############################################################################
## package 'secr'
## rbind.mask.R
## 2017-07-25 moved from methods.R
#############################################################################

rbind.mask <- function (...) {
    # combine 2 or more mask objects
    
    ##    no check for multi-session masks at present
    ##        stop ('rbind of multi-session mask not implemented')
    ##
    dropduplicates <- TRUE   ## always
    allargs <- list(...)
    spacing <- attr(allargs[[1]],'spacing')
    area    <- attr(allargs[[1]], 'area')
    check <- function (x) {
        if (!is(x,'mask'))
            stop ("arguments must be mask objects")
        if (attr(x,'spacing') != spacing)
            stop ("arguments must have same 'spacing' attribute")
        if (attr(x,'area') != area)
            stop ("arguments must have same area attribute")
    }
    sapply (allargs, check)
    temp <- rbind.data.frame(...)
    class(temp) <- c('mask', 'data.frame')
    tempcov <- lapply(allargs, covariates)
    covariates(temp) <- do.call(rbind, tempcov)  ## pass list of dataframes
    
    if (dropduplicates) {
        dupl <- duplicated(temp)
        droppedrows <- sum(dupl)
        if (droppedrows>0) {
            covariates(temp) <- covariates(temp)[!dupl,]
            temp <- temp[!dupl,]
            warning (droppedrows, " duplicate points dropped from mask")
        }
    }
    
    attr(temp,'type')        <- 'rbind'
    attr(temp,'meanSD')      <- getMeanSD(temp)
    attr(temp,'area')        <- area
    attr(temp,'spacing')     <- spacing
    xl <- range(temp$x) + spacing/2 * c(-1,1)
    yl <- range(temp$y) + spacing/2 * c(-1,1)
    
    attr(temp,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    temp
}
############################################################################################
