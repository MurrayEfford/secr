############################################################################################
## package 'secr'
## model.average.R
## 2023-09-30 defunct
############################################################################################

model.average <- function (..., realnames = NULL, betanames = NULL,
    newdata = NULL, alpha = 0.05, dmax = 10, covar = FALSE, average =
        c('link', 'real'), criterion = c('AICc','AIC'), CImethod =
        c('Wald', 'MATA')) {
    
    # .Deprecated("modelAverage", package="secr", "This function has been renamed.", 
    #     old = as.character(sys.call(sys.parent()))[1L])

    # 2023-09-30    
    .Defunct("modelAverage", package="secr", "This function has been renamed.")
}
############################################################################################

