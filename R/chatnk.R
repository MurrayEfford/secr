###############################################################################
## package 'secr'
## chatnk.R
## 2022-11-18, 29
## 2023-04-26
###############################################################################

Fletcher.chat <- function(observed, expected, np) {
    K <- length(observed)
    X2 <- sum((observed - expected)^2 / expected)
    si <- sum((observed - expected) / expected) / K
    nu <- K-np
    list(
        expected = expected, 
        observed = observed, 
        stats = c(
            mean.expected = mean(expected), 
            var.expected = sd(expected)^2,
            mean.observed = mean(observed), 
            var.observed = sd(observed)^2, 
            si = si,
            nu = nu,
            cX2 = X2/nu),
        chat = X2/nu / (1 + si)
    )
}

chat.nk.sess <- function(object, D, capthist, mask, detpar) {
    
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
    
    np <- length(object$betanames)
    if (np > (nrow(traps(capthist))-1)) stop ("c-hat not estimated when np > K-1")
    
    Fletcher.chat(nk, expected.nk, np)
    
}

chat.nk <- function(object) {
    det <- unlist(detector(traps(object$capthist)))
    if (!all(det %in% c('multi','proximity','count'))) {
        stop("chat.nk available only for multi, proximity and count detectors")
    }
    
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
        getdet <- function(x) {
            ok <- rownames(x) %in% c('g0','lambda0','sigma','z')
            as.list(setNames(x[ok,'estimate'], rownames(x)[ok]))
        }
        detparlist <- lapply(predict(object), getdet)
        mapply(chat.nk.sess, 
            D = Dlist, 
            capthist = object$capthist,   # 2023-04-23 previously object$capthist[[1]] 
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
            pred <- predict(object)
            D <- pred['D', 'estimate']
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

# 2022-11-20
# experimental adjustment of SE and CL
# apply to density linear predictor (on link scale)
# include in predict.secr?

adjustVarD <- function(object, chatmin = 1.0, alpha = 0.05, chat = NULL) {
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
    if (!inherits(object, 'secr')) {
        # expect dataframe input
        pred <- list(object['D',])  
        if (is.null(chat)) stop ("specify chat for data.frame input")
    }
    else {
        if ('D' %in% names(object$model)) {
            pred <- predict(object)
        }
        else {
            pred <- derived(object, se.esa = FALSE, se.D = TRUE)
        }
        
        # pred and chat should be lists with one component per session
        if (ms(object)) {
            if (is.null(chat)) chat <- sapply(chat.nk(object), '[[', 'chat') # vector of chat
            pred <- lapply(pred, '[', 'D',)
        }
        else {
            if (is.null(chat)) chat <- chat.nk(object)$chat
            pred <- list(pred['D',])
        }
    }
    do.call(rbind, mapply(adjustonesession, pred, chat, 
        MoreArgs = list(chatmin = chatmin, alpha = alpha),
        SIMPLIFY = FALSE))
}