# 2022-11-20, 2025-05-20

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
        
        # chat should be vector with one element per session
        if (is.null(chat)) {
            # chat <- unlist(chat.nk(object, verbose = FALSE, type = 'Fletcher')) 
            stop ("secr >= 5.2.2 adjustVarD() requires chat to be provided; nk default no longer offered")
        }
        # pred should be list with one component per session
        if (ms(object)) {
            pred <- lapply(pred, '[', 'D',)
        }
        else {
            pred <- list(pred['D',])
        }
    }
    # if (length(pred) != length(chat)) stop ("mismatch of chat vector and predicted values")
    do.call(rbind, mapply(adjustonesession, pred, chat, 
                          MoreArgs = list(chatmin = chatmin, alpha = alpha),
                          SIMPLIFY = FALSE))
}
