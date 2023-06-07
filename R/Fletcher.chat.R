###############################################################################
## package 'secr'
## Fletcher.chat.R
## 2023-05-05
## 2023-05-10 switch order Fletcher, Wedderburn when type = 'both'
###############################################################################

# Fletcher.chat is called by chat.nk.sess;

Fletcher.chat <- function(observed, expected, np, verbose = TRUE, 
    type = c('Fletcher', 'Wedderburn', 'both'), multinomial = FALSE) {
    type <- match.arg(tolower(type), choices = c('fletcher', 'wedderburn', 'both'))
    if (is.list(observed)) {
        # apply Fletcher.chat recursively to each component of 'observed'
        if (type == 'both') {
            list(
                Fletcher = sapply(observed, Fletcher.chat, expected = expected, 
                    np = np, verbose = FALSE, type = 'Fletcher'),
                Wedderburn = sapply(observed, Fletcher.chat, expected = expected, 
                    np = np, verbose = FALSE, type = 'Wedderburn')
            )
        }
        else {
            sapply(observed, Fletcher.chat, expected = expected, 
                np = np, verbose = FALSE, type = type)
        }
    }
    else {
        K <- length(observed)
        nu <- K - np - multinomial
        protected <- TRUE
        
        # Wedderburn
        cX2 <- sum( (observed - expected)^2 / expected ) / nu  
        if (protected) {
            # Guard against divide-by-zero 2023-06-07
            ok <- observed>0
            if (any(observed>0 & expected<1e-20)) {
                warning ("observed>0 when expected<1e-20")
                cX2 <- NA
            }
            else {
                cX2 <- (
                    sum(observed[ok]^2 / expected[ok])
                    - 2 * sum(observed)
                    + sum(expected)
                ) / nu
            }
        }
        # Fletcher
        sbar <- sum( (observed - expected) / expected ) / (K - multinomial)
        if (protected) {
            sbar <- (sum( observed[ok] / expected[ok] ) ) / (K - multinomial) - 1
        }
        chat <- cX2 / (1 + sbar)  
        if (verbose) {
            list(
                expected = expected, 
                observed = observed, 
                stats = c(
                    mean.expected = mean(expected), 
                    var.expected = sd(expected)^2,
                    mean.observed = mean(observed), 
                    var.observed = sd(observed)^2, 
                    sbar = sbar,
                    nu = nu,
                    cX2 = cX2),
                chat = chat
            )
        }
        else {   # scalar
            if (type == "fletcher") chat 
            else if (type == "wedderburn") cX2
            else c(Wedderburn = cX2, Fletcher = chat) # "both"
        }
    }
}
