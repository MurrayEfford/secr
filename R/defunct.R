############################################################################################
## package 'secr'
## defunct.R
## 2023-10-18
############################################################################################

model.average <- function (..., realnames = NULL, betanames = NULL,
                           newdata = NULL, alpha = 0.05, dmax = 10, covar = FALSE, average =
                               c('link', 'real'), criterion = c('AICc','AIC'), CImethod =
                               c('Wald', 'MATA')) {
    .Defunct("modelAverage", package="secr", "This function has been renamed modelAverage.")
}
############################################################################################

ip.secr <- function (...) {
    .Defunct("", package="secr", "Function ip.secr has been removed. See ipsecr.fit() in package 'ipsecr'.")
}
############################################################################################
# 
# pfn <- function (...) {
#     .Defunct("", package="secr", "Function pfn has been removed. See package 'ipsecr'.")
# }
# ############################################################################################
# make.newdata <- function (...) {
#     .Defunct("makeNewData", package="secr", "This function has been renamed makeNewData.")
# }
# ############################################################################################
# secr.make.newdata <- function (...) {
#     .Defunct("makeNewData", package="secr", "This function has been renamed makeNewData.")
# }
# ############################################################################################
# read.SPACECAP <- function (...) {
#     .Defunct("", package="secr", "Function read.SPACECAP has been removed.")
# }
# ############################################################################################
# write.SPACECAP <- function (...) {
#     .Defunct("", package="secr", "Function write.SPACECAP has been removed.")
# }
# ############################################################################################
