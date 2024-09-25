#############################################################################
## package 'secr'
## xy2CH.R
## 2024-09-25 moved from utility.R as exists only for export
#############################################################################

## convert telemetryxy attribute of a combined dataset into a standalone capthist

## TO DO: option of telemetry or polygon output

xy2CH <- function (CH, inflation = 1e-8) {
    if (ms(CH)) {
        out <- lapply(CH, xy2CH, inflation)
        class(out) <- c('capthist', 'list')
        out
    }
    else {
        xylist <- telemetryxy(CH)
        if (is.null(xylist))
            stop ("requires 'telemetryxy' attribute")
        n <- length(xylist)
        neach <- sapply(xylist, nrow)
        allxy <- do.call(rbind, xylist)
        
        trps <-  allxy[chull(allxy),]
        trps <- rbind(trps, trps[1,,drop=F])
        trps <- inflate(trps, 1 + inflation)  ## see also telemetry.R
        
        trps <- as.data.frame(trps)
        dimnames(trps) <- list(1:nrow(trps), c('x','y'))
        class(trps) <- c("traps","data.frame")
        detector(trps) <- "polygon"
        polyID(trps) <- factor(rep(1,nrow(trps)))
        
        rown <- rep(names(xylist), neach)
        newCH <- array(neach, dim = c(n, 1, 1))
        attr(newCH, "detectedXY") <- allxy
        if (!is.null(covariates(CH))) {
            rowlookup <- match(names(xylist), rownames(CH))
            covariates(newCH) <- covariates(CH)[rowlookup,, drop=FALSE]
        }
        class(newCH) <- "capthist"
        traps(newCH) <- trps
        newCH
    }
}

