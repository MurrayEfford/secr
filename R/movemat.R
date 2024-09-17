movemat <- function (CH) {
    if (ms(CH)) { lapply(CH, movemat)}
    else {
        nocc <- ncol(CH)
        ntraps <- nrow(traps(CH))
        if (nocc<2) stop ("movemat requires multiple occasions")
        df <- as.data.frame(CH)
        df$trapno <- as.numeric(factor(df$TrapID))
        byID <- split(df, df$ID)
        oneID <- function (df1) {
            movematcpp(ntraps, as.integer(df1$trapno))
        }
        tmp <- lapply(byID, oneID)
        tmp
    }
}
# secr:::movemat(captdata)