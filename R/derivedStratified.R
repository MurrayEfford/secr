derivedStratified <- function (object, strata, mask = NULL, sessnum = NULL, plt = FALSE, ...) {

    if (ms(object)) {
        ## recursive call if MS
        stop("derivedStratified not yet implemented for multi-session fit")
        sessnames <- session(object$capthist)
        nsess <- length(sessnames)
        output <- vector('list', nsess )
        for (i in 1:nsess) {
            output[[i]] <- derivedStratified(object, strata, sessnum = i, mask, ...)
        }
        names(output) <- sessnames
        output
    }
    else {

        if (is.null(sessnum)) {
            capthist <- object$capthist
            sessnum <- 1
        }
        else
            capthist <- object$capthist[[sessnum]]

        if (is.null(mask))  ## multisession case?
            mask <- object$mask

        CH <- reduce(capthist, outputdetector = 'proximity')
        trapstrata <- group.factor(traps(CH), groups=strata)

        ## step 1: get stratum-specific counts of individuals nh
        nh <- tapply(apply(abs(CH), 3, sum), trapstrata, sum)
        if (nrow(CH) != sum(nh))
            warning ("some animals detected in more than one stratum")

        ## step 2: get stratum-specific a and D-hat
        ## WARNING: HAVE NOT ALLOWED
        ## INDIVIDUAL-SPECIFIC OR TRAP-SPECIFIC MODELS
        ## IS SUBSETTING PIA REQUIRED?
        stratalevels <- levels(trapstrata)
        nstrata <- length(stratalevels)
        if (is.null(mask))
             mask <- object$mask

        which.trap <- nearesttrap(mask, traps(CH))
        out <- vector('list', nstrata+1)
        for (i in 1:nstrata) {
            stratumi <- stratalevels[i]
            object$capthist <- subset(CH, traps = trapstrata == stratumi)
            if (all (strata %in% names(covariates(object$mask)))) {
                ## derive mask from covariates
                maskstrata <- group.factor(mask, groups = strata) ## not yet
            }
            else {
                maskstrata <- trapstrata[which.trap]
            }

            object$mask <- subset(mask, maskstrata == stratumi)
            tmp <- derived(object, se.esa = TRUE, se.D = TRUE, ...)
            out[[i]] <- data.frame(esa = tmp[1,1], SE.esa = tmp[1,2],
                                   wh = nrow(object$mask), tmp[2,])
        }
        out <- do.call(rbind, out)
        out$wh <- out$wh/sum(out$wh)
        out2 <- rbind(out, rep(NA,ncol(out)))
        rownames(out2) <- c(stratalevels,'Total')
        out2[nstrata+1, 'wh'] <- sum(out$wh)
        out2[nstrata+1, 'estimate'] <- sum(out$wh * out$estimate)
        varD <- sum(out$wh^2 * out$SE.estimate^2) ## IGNORING ANY COVARIANCES
        out2[nstrata+1, 'SE.estimate'] <- varD^0.5
        out2[nstrata+1, 'CVD'] <- out2[nstrata+1, 'SE.estimate'] / out2[nstrata+1, 'estimate']

        if (plt) {
            covariates(mask) <- data.frame(stratum = )
        }

        out2

    }
}
####################################################################################

