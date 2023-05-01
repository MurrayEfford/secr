###############################################################################
## package 'secr'
## chatnk.R
## 2022-11-18, 29
## 2023-04-26 to 2023-05-01
###############################################################################

# Fletcher.chat is called by chat.nk.sess;
# it is not exported by secr
Fletcher.chat <- function(observed, expected, np, verbose = TRUE, 
    type = c('Fletcher', 'Pearson', 'both')) {
    type <- match.arg(tolower(type), choices = c('fletcher', 'pearson', 'both'))
    if (is.list(observed)) {
        # apply Fletcher.chat recursively to each component of 'observed'
        if (type == 'both') {
            list(
                Fletcher = sapply(observed, Fletcher.chat, expected = expected, 
                    np = np, verbose = FALSE, type = 'Fletcher'),
                Pearson = sapply(observed, Fletcher.chat, expected = expected, 
                    np = np, verbose = FALSE, type = 'Pearson')
            )
        }
        else {
            sapply(observed, Fletcher.chat, expected = expected, 
                np = np, verbose = FALSE, type = type)
        }
    }
    else {
        K <- length(observed)
        X2 <- sum((observed - expected)^2 / expected)
        si <- sum((observed - expected) / expected) / K
        nu <- K-np
        cX2 <- X2 / nu
        chat <- cX2 / (1 + si)
        if (verbose) {
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
                    cX2 = cX2),
                chat = chat
            )
        }
        else {   # scalar
            if (type == "fletcher") chat 
            else if (type == "pearson") cX2
            else c(Fletcher = chat, Pearson = cX2) # "both"
        }
    }
}

chat.nk.sess <- function(object, D, capthist, mask, detpar, nsim, ...) {
    
    ## c-hat for one session
    
    ## potential development:
    ## check PIA0 for variation over occasions, animals, detectors?
    ## substitute detectpar by lookup of Xrealparval0 with PIA0?
    
    noccasions <- dim(capthist)[2]
    traps <- traps(capthist)
    expected.nk <- Enk(
        D          = D, 
        mask       = mask, 
        traps      = traps, 
        detectfn   = object$detectfn, 
        detectpar  = detpar, 
        noccasions = noccasions,
        binomN     = object$binomN, 
        userdist   = object$details$userdist, 
        ncores     = NULL,
        nrepl      = NULL)   # do not simulate expected
    
    np <- length(object$betanames)
    if (np > (nrow(traps)-1)) stop ("c-hat not estimated when np > K-1")
    
    if (is.null(nsim) || nsim < 1) {
        observed.nk <- apply(apply(abs(capthist),c(1,3),sum)>0, 2, sum)
        
    }   
    else {
        # simulate a list of 'observed' nk vectors
        onesimnk <- function (r) {
            pop <- sim.popn(D, core = mask, model2D = 'IHP')
            ch <- sim.capthist(
                traps      = traps, 
                popn       = pop, 
                detectfn   = object$detectfn, 
                detectpar  = detpar, 
                noccasions = noccasions, 
                nsessions  = 1, 
                binomN     = object$binomN)
            apply(apply(ch, c(1,3), max),2,sum)  # individuals per detector
        }
        observed.nk <- lapply(1:nsim, onesimnk)
    }
    
    Fletcher.chat(observed.nk, expected.nk, np, ...)
    
}

chat.nk <- function(object, nsim = NULL, ...) {
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
            MoreArgs = list(object = object, nsim = nsim, ...), 
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
        chat.nk.sess (object, D, capthist, mask, detpar, nsim = nsim, ...)
        
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