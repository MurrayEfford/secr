###############################################################################
## package 'secr'
## onAttach.R
## last changed 2011-06-16 2013-04-20 2016-10-07 2024-09-28 2024-10-02
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- paste0(packageVersion('secr'), .localstuff$packageType)
    packageStartupMessage( 
        "This is secr ", version, ". For overview type ?secr"
        )
}

## .onLoad is preferred if actions are required for single functions 
## that may be called without attaching package
