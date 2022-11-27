###############################################################################
## package 'secr'
## chatnk.R
## 2022-11-18
###############################################################################

chat.nk.sess <- function(object, D, capthist, mask, detpar) {
    
    if (ms(object)) {
        stop("expect single session")
    }
    else {
        ## check PIA0 for variation over occasions, animals, detectors?
        ## substitute detectpar by lookup of Xrealparval0 with PIA0?
        
        expected.nk <- Enk(
            D          = D, 
            mask       = mask, 
            traps      = traps(capthist), 
            detectfn   = object$detectfn, 
            detectpar  = detpar, 
            noccasions = dim(capthist)[2],
            binomN     = object$binomN, 
            userdist   = object$details$userdist, 
            ncores     = NULL)
        
        nk <- apply(apply(abs(capthist),c(1,3),sum)>0, 2, sum)
        X2 <- sum((nk - expected.nk)^2 / expected.nk)
        si <- mean((nk - expected.nk) / expected.nk)
        nu <- nrow(traps(capthist)) - 1
        chat <- X2/nu / (1 + si)
        list(
            expected.nk = expected.nk, 
            nk = nk, 
            stats = c(
                expected.nk = mean(expected.nk), 
                var.expected = sd(expected.nk)^2,
                mean.nk = mean(nk), 
                var.nk = sd(nk)^2, 
                cX2 = X2/nu),
            chat = chat
        )
        ##-------------------------------------------------------------
        
    }
}

chat.nk <- function(object) {
    
    if (ms(object)) {
        if (object$CL) {
            Dlist <- lapply(derived(object), '[', 'D', 'estimate')
        }
        else if (object$model$D == ~1) {
            Dlist <- lapply(predict(object), '[', 'D', 'estimate')
        }
        else {
            Dlist <- lapply(covariates(predictDsurface(object)), '[[', 'D.0')
        }
        detparlist <- detectpar(object)
        ch <- object$capthist
        object$capthist <- object$capthist[[1]] # to fool not ms
        mapply(chat.nk.sess, 
            D = Dlist, 
            capthist = ch, 
            mask = object$mask, 
            detpar = detparlist, 
            MoreArgs = list(object = object), 
            SIMPLIFY = FALSE)
        
    }
    else {
        
        ##-------------------------------------------------------------
        ## Restrict application
        if (!inherits(object, 'secr'))
            stop ("chat.nk expects fitted secr model")
        if (length(object$fixed)>0)
            stop ("chat.nk does not yet accept models with one or more fixed parameters")
        if (!is.null(object$groups))
            stop ("chat.nk does not yet accept models with grouping")
        if (object$details$nmix>1)
            stop ("chat.nk does not yet handle mixtures")
        if (length(table(object$design0$PIA))>1 || length(table(object$design$PIA))>1)
            stop ("chat.nk does not yet handle varying detection probabilities")
        
        ##-------------------------------------------------------------
        
        ## Density (scalar or length = nrow(object$mask))
        
        if (object$CL) {
            D <- derived(object)['D', 'estimate']
        }
        else if (object$model$D == ~1) {
            D <- predict(object)['D', 'estimate']
        }
        else {
            D <- covariates(predictDsurface(object))$D.0
        }
        ##-------------------------------------------------------------
        
        capthist <- object$capthist
        mask <- object$mask
        detpar <- detectpar(object)
        chat.nk.sess (object, D, capthist, mask, detpar)
        
    }
}

# sapply(chat.nk(oven), '[[','chat')

# 2022-11-20
# unpublished
# experimental adjustment of SE and CL
# apply to density linear predictor (on link scale)

# also derived.secr?
# include in predict.secr?

chatVariance <- function(object, chatmin = 1.0, alpha = 0.05) {
    adjustonesession <- function (pred, chat, chatmin = 1.0, alpha = 0.05) {
        link <- pred['D', 'link']
        if (is.null(link)) link <- 'log'
        D    <- pred['D', 'estimate']
        seD  <- pred['D', 'SE.estimate']
        selinkD <- se.transform(D, seD, link) 
        linkD   <- transform(D, link)
        chat <- max(chatmin, chat)
        z <- abs(qnorm(1-alpha/2))  
        pred['D', c('lcl','ucl')] <- untransform(linkD + z*selinkD*c(-1,1)*chat^0.5, link)
        pred['D', c('lcl','ucl')] <- untransform(linkD + z*selinkD*c(-1,1)*chat^0.5, link)
        pred['D', 'SE.estimate']  <- se.untransform(linkD, selinkD*chat^0.5, link)
        pred['D', 'c-hat'] <- chat
        pred
    }
    if ('D' %in% names(object$model)) {
        pred <- predict(object)
    }
    else {
        pred <- derived(object, se.esa = FALSE, se.D = TRUE)
    }
    
    # pred and chat should be lists with one component per session
    if (ms(object)) {
        pred <- lapply(pred, '[', 'D',)
        chat <- sapply(chat.nk(object), '[[', 'chat') # vector of chat
    }
    else {
        pred <- list(pred['D',])
        chat <- chat.nk(object)$chat
    }
    do.call(rbind, mapply(adjustonesession, pred, chat, 
        MoreArgs = list(chatmin = chatmin, alpha = alpha),
        SIMPLIFY = FALSE))
}

#-------------------------------------------------------------------------------
# cf Paul Hiemstra Aug 23 2012 StackOverFlow
binarystring <- function(number, noBits) {
    binary_vector = rev(as.numeric(intToBits(number)))
    paste(binary_vector[-(1:(length(binary_vector) - noBits))], collapse='')
}

# steps towards overdispersion of detection Nov 2022

chat.i <- function (object) {
    ch <- object$capthist
    noccasions <- ncol(ch)
    chis <- apply(abs(ch),1:2,sum)
    chis <- pmin(chis,1) # ignore detection at multiple traps
    observedCH <- apply(chis,1,paste, collapse='')
    possibleCH <- sapply(1:2^noccasions, binarystring,noccasions)
    table(factor(observedCH, levels = possibleCH))  
}
# chat.i(secrdemo.0)
#-------------------------------------------------------------------------------
