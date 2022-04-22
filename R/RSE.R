## extract relative standard error from fitted secr or openCR model
## default is vector for all parameters in model

## 2020-05-06 RSE may return matrix of values

RSE <- function (fit, parm = NULL, newdata = NULL) {
    if (is.null(parm)) {
        parm <- names(fit$parindx)
    }
    if (length(parm) > 1) {
        sapply(parm, RSE, fit = fit, newdata = newdata)
    }
    else {
        ## condition ensures 
        ## (i) parm in model (ii) no modelled variation in parm
        pindex <- length(fit$parindx[[parm]])
        if (pindex == 0) {
            warning ("parameter ", parm, " not in model")
        }
        else if (pindex == 1 && fit$link[[parm]] == 'log') {
            if (!is.null(newdata)) {
                warning ("non-null newdata ignored for ", parm)
            }
            ## Efford & Boulanger MEE 2019
            sqrt(exp(vcov(fit)[parm, parm])-1)    
        }
        else {
            if (pindex > 1 && is.null(newdata)) {
                warning ("parameter ", parm, " varies in model; consider specifying newdata")
            }
            pred <- predict(fit, newdata = newdata)
            if (parm %in% rownames(pred)) {
                pred[parm, 'SE.estimate'] / pred[parm, 'estimate']
            }
            else if (parm %in% rownames(pred[[1]])) {
                # modified 2020-05-06 to return all
                getone <- function (x) {
                    x[parm, 'SE.estimate'] / x[parm, 'estimate']
                }
                sapply(pred, getone)
            }
            else if (parm %in% names(pred)) {
                pred[[parm]][1,'SE.estimate'] / pred[[parm]][1,'estimate']
            }
            else NA
        }
    }
}