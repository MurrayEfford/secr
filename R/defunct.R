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

# 2024-09-19
# derivedSystematic( object, xy, design = list(), basenx = 10, df = 9, extrapolate = TRUE,
#                    alpha = 0.05, loginterval = TRUE, independent.esa = FALSE, keep = FALSE,
#                    ncores = NULL)

# \item{design}{list specifying systematic design (see Details)}
# \item{basenx}{integer number of basis grid points in x-dimension}
# \item{df}{integer number of degrees of freedom for gam}
# \item{extrapolate}{ logical; if FALSE then boxlet p values are inferred from nearest point inside convex hull of grid}
# \item{keep}{logical; if TRUE then derivedSystematic saves key intermediate values as attributes}
# \item{ncores}{integer}
 
# \code{derivedSystematic} implements the 'boxlet' variance estimator of Fewster (2011) for systematic designs using clustered detectors (an alternative to \code{derivedCluster} and \code{derivedSessions}). The method is experimental in secr 3.2.0 and may change. The `design' argument is a list with components corresponding to arguments of \code{\link{make.systematic}}, (\code{n} and \code{origin} are ignored if provided):
# 
# \tabular{llll}{
# Component \tab Description \cr
# \code{cluster} \tab traps object for a single cluster \cr
# \code{region} \tab 2-column matrix or SpatialPolygons \cr
# \code{spacing} \tab spacing between cluster origins \cr
# \code{...} \tab other arguments passed to \code{\link{trap.builder}} \cr
# \tab e.g. \code{edgemethod}, \code{exclude}, \code{exclmethod} \cr
# }
# 
# If \code{region} is omitted from \code{design} then an attempt will be made to retrieve it from the mask attribute of \code{object} (this works if the call to \code{\link{make.mask}} used \code{keep.poly = TRUE}).

# ############################################################################################

