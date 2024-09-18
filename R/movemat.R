movemat <- function (CH, randomize = TRUE) {
    if (ms(CH)) { lapply(CH, movemat)}
    else {
        n <- nrow(CH)
        nocc <- ncol(CH)
        ntraps <- nrow(traps(CH))
        if (nocc<2) stop ("movemat requires multiple occasions")
        df <- as.data.frame(CH)
        # ensure trap order does not change
        df$trapno <- as.numeric(factor(df$TrapID, levels = rownames(traps(CH))))
        byID <- split(df, df$ID)
        oneID <- function (df1) {
            k <- df1$trapno
            if (randomize) {
                k <- sample(k, length(k), replace = FALSE)
            }
            movematcpp(as.integer(ntraps), as.integer(k))
        }
        tmp <- lapply(byID, oneID)
        array(unlist(tmp), dim = c(ntraps, ntraps, n))
    }
}

# secr:::movemat(captdata)