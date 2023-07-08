############################################################################################
## package 'secr'
## LRtest.R
## likelihood ratio test for two models (assumed to be nested)
## last changed 2011 12 13 to include openCR
## 2013-10-29 fixed name of p.value; use logLik
## 2017-10-14 allow for length class vector > 1 [any(class...)]
## 2017-10-14 generalised by getting np from logLik
## 2023-06-26 class(model1)[1]
############################################################################################

LR.test <- function (model1, model2) {
    if (is.null(getS3method('logLik', class(model1)[1], TRUE)))
        stop ("no logLik method for model")
    if (any(class(model1) != class(model2)))
        stop ("models must have same class")
    ## 2019-11-29
    if (inherits(model1, "secr") && !all(AICcompatible(model1, model2))) {
        stop ("models not compatible for LR.test")
    }
    
    statistic <- as.numeric(2 * abs(logLik(model1) - logLik(model2)))
    if (length(statistic) != 1)
        stop ("problem with 'model1' or 'model2'")
    
    parameter <- abs(attr(logLik(model1), "df") -  attr(logLik(model2), "df"))
    p.value <- 1 - pchisq(statistic, parameter)
    names(statistic) <- 'X-square'
    names(parameter) <- 'df'
    names(p.value) <- NULL
    s1name <- deparse(substitute(model1))
    s2name <- deparse(substitute(model2))
    temp <- list (
        statistic = statistic,
        parameter = parameter,
        p.value = p.value,
        method = 'Likelihood ratio test for two models',
        data.name = paste (s1name, 'vs', s2name)
    )
    class(temp) <- 'htest'
    temp
}
############################################################################################

