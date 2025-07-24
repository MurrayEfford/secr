## summary.secr.R
## 2018-10-14

summary.secr <- function (object, newdata = NULL, alpha = 0.05, deriv = FALSE, ...) {

    out <- vector('list')

    # if (call) {
    #     out$call <- object$call
    # }

    out$versiontime <- paste0(object$version, ', run ', object$starttime, ', elapsed ', round(object$proctime,2), ' s')
    if (!is.null(object$details$newdetector)) {
        out$newdetector <- object$details$newdetector
    }

    ###################
    ## Data description

    if (ms(object$capthist)) {
        out$capthist <- summary(object$capthist, terse = TRUE, moves = TRUE, tpa = TRUE)
        if (!is.null(markocc(traps(object$capthist)[[1]]))) {
            out$Tu <- sapply(object$capthist, function(y) attr(y,'Tu',exact = TRUE))
        }
        out$trapsummary <- summary(traps(object$capthist))
        out$detector <- sapply(detector(traps(object$capthist)), '[',1)
        ## redundant because of chsummary
        out$xyl <- NULL
    }
    else {
        trp <- traps(object$capthist)
        out$traps <- data.frame (Detector = detector(trp)[1],
                                 Number = nrow(trp),
                                 Spacing = spacing(trp))
        if (is.null(usage(trp)))
            out$traps$UsagePct <- 100
        else
            out$traps$UsagePct <- 100 * sum(usage(trp))/length(usage(trp))
        if (length(detector(trp))>1)
            out$detector <- detector(trp)
        out$capthist <- summary(object$capthist, terse = TRUE, moves = TRUE, tpa = TRUE)

        if (!is.null(markocc(traps(object$capthist) ))) {
            defaultmarkresight <- list(Tu='as.is', Tm='as.is', Tn='ignore')
            markresightcontrol <- secr_replacedefaults(defaultmarkresight, object$details$markresight)
            out$Tu <- Tu(object$capthist)
            out$Tm <- Tm(object$capthist)
            out$Tn <- Tn(object$capthist)
            out$unresolvedocc <- markocc(traps(object$capthist)) == -1
            out$unmarkedocc <- markocc(traps(object$capthist)) == 0

            # if (!(is.null(Tu) | markresightcontrol$Tu=='ignore'))
            #     cat ('N unmarked sght :  ', sum(unlist(Tu)), ' (c-hat ', round(object$details$chat[1],3),  ')\n', sep='')
            # if (!(is.null(Tm) | markresightcontrol$Tm=='ignore'))
            #     cat ('N nonID sghting :  ', sum(unlist(Tm)), ' (c-hat ', round(object$details$chat[2],3), ')\n', sep='')
            # if (!(is.null(Tn) | markresightcontrol$Tn=='ignore'))
            #     cat ('N other sghting : ', sum(unlist(Tn)), '\n')
        }

        if ('g' %in% object$vars) {
            Groups  <- table(secr_group.factor(object$capthist, object$groups))
            temp <- paste (names(Groups), Groups, collapse=', ', sep='=')
            out$groupinfo <- paste('(',temp,')', sep='')
        }

        # out$noccasions <- ncol(object$capthist)

        ## Telemetry
        # xyl <- telemetryxy(object$capthist)
        # if (!is.null(xyl)) {
        #     ## zeros <- sum(apply(abs(object)>0,1,sum)==0)
        #     ntelem <- sapply(xyl, nrow)
        #     nteldet <- if ((nrow(object$capthist) == 0) | (all(detector(traps(object$capthist))=='telemetry')))
        #         0 else
        #             sum(apply(abs(object$capthist)[,-ncol(object$capthist),,drop=FALSE]>0,1,any) [row.names(object$capthist) %in% names(xyl)])
        #     ## cat ('Known all-zero  : ', zeros, '\n')
        #     cat ('Telemetry       : ', length(xyl), 'animals,', nteldet, 'detected\n')
        #     cat ('Telemetry locns : ', paste(range(ntelem), collapse='-'), 'per animal (mean',
        #          round(mean(ntelem),2), 'sd', paste(round(sd(ntelem),2), ')',sep=''), '\n')
        # }


    }   # end of single-session
    if (any(out$detector %in% .localstuff$countdetectors)) {
        out$countmodel <- if (object$details$binomN == 0) 'Poisson'
        else if (object$details$binomN == 1) 'Binomial, size from usage'
        else if (object$details$binomN < 0) paste0('Negative binomial k = ', abs(object$details$binomN))
        else if (object$details$binomN > 1) paste0('Binomial', object$details$binomN)
    }

    if (ms(object$mask)) {
        out$mask <- data.frame(Cells = sapply(object$mask, nrow),
                                      Spacing = sapply(object$mask, spacing))
        if (length(maskarea(object$mask[[1]]))==0)
            out$mask$Length <- sapply(object$mask, masklength)
        else
            out$mask$Area <- sapply(object$mask, maskarea)
    }
    else {
        out$mask <- data.frame(Cells = nrow(object$mask), Spacing = spacing(object$mask))
        if (length(maskarea(object$mask))==0)
            out$mask <- cbind(out$mask, Length = masklength(object$mask))
        else
            out$mask <- cbind(out$mask, Area = maskarea(object$mask))
    }

    ####################
    ## Model description

    out$modeldetails <- data.frame(CL = as.character(object$CL),
                                   fixed = fixed.string(object$fixed),
                                   distribution = if (!object$CL) object$details$distribution else '',
                                   hcov = if (!is.null(object$hcov)) object$hcov else '',
                                   relativeD = as.character(object$details$relativeD))

    out$AICtable <- AIC(object)

    if (!is.null(object$details$userdist)) {
        out$userdist <- if (is.matrix(object$details$userdist))
             'static (matrix)'
        else if (is.function(object$details$userdist))
            'dynamic (function)'
        else ''
    }

    out$coef <- coef(object)
    # out$beta.vcv <- object$beta.vcv

    # # scale newdata covariates... NOT FINISHED 10 05 08
    # meanSD <- attr(object$mask,'meanSD',exact = TRUE)
    # if (!is.null(newdata)) {
    #     for (i in 1:length(newdata)) {
    #         ind <- match (names(newdata[i]),names(meanSD))
    #         if (ind>0 & !is.na(meanSD[1,ind]))
    #             newdata[[i]] <- (newdata[[i]] - meanSD[1,ind]) / meanSD[2,ind]
    #     }
    # }

    out$predicted <- predict (object, newdata, type = "response", alpha = alpha)

    #################################
    # Derived parameters
    #################################
    if (deriv) {
        out$derived <- derived(object, alpha=alpha, ...)
    }

    ## remove distracting row names
    for (i in 1:length(out)) {
        if (is.data.frame(out[[i]]))
            if (nrow(out[[i]])==1 & (!names(out)[i] %in% c('coef', 'derived', 'predicted')))
                rownames (out[[i]]) <- ''
    }
    class(out) <- c("summary.secr")
    out
}
############################################################################################

print.summary.secr <- function (x, ...) {
    class(x) <- NULL
    print(x)
}
############################################################################################

AIC.summary.secr <- function (object, ..., sort = TRUE, k = 2, dmax = 10, 
                              criterion = c('AIC','AICc')) {
    criterion <- match.arg(criterion)
    allargs <- list(object, ...)
    output <- do.call(rbind, lapply(allargs, '[[', "AICtable"))
    rownames(output) <- NULL
    output$delta <- output[,criterion] - min(output[, criterion])
    ## optional sort
    if (sort) output <- output[order(output$delta),]
    ## AICwt with dmax
    OK <- abs(output$delta) < abs(dmax)
    sumdelta <- sum(exp(-output$delta[OK]/2))
    output$wt <- ifelse ( OK, round(exp(-output$delta/2) / sumdelta,4), 0)
    names(output)[7] <- paste('d',criterion,sep='')
    names(output)[8] <- paste(criterion,'wt',sep='')
    if (nrow(output)==1) { output[,8] <- NULL; output[,7] <- NULL}
    output
}
############################################################################################
