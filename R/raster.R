## register S3 classes as (virtual) S4 classes
setOldClass('mask')
setOldClass('Dsurface')

## 'raster' and 'terra' are not imported, because loading them is slow and secr does not
## otherwise need them. The 'raster' and 'rast' methods for mask and Dsurface are therefore
## registered when (and if) those packages are loaded; see .onLoad. The methods are stored in a
## private environment because a loaded namespace is locked.
.rastermethods <- new.env()

## define 'raster' methods for mask and Dsurface
registerRasterMethods <- function () {
    methods::setMethod(raster::raster, methods::signature(x = "mask"),
        ## define internal function for coercing mask to raster
        function(x,  covariate, values = 1, crs = NA) {
            if (missing(covariate)) {
                covariate <- 'tempcov'
                covariates(x) <- data.frame(tempcov = rep(values, length.out=nrow(x)))
            }
            mask <- rectangularMask(x)
            bbox <- attr(x, 'boundingbox')
            xmx <- max(bbox$x)
            xmn <- min(bbox$x)
            ymx <- max(bbox$y)
            ymn <- min(bbox$y)
            nx <- round((xmx - xmn) / spacing(mask))
            ny <- round((ymx - ymn) / spacing(mask))
            tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
            mat <- matrix (tmp, nrow = ny, ncol = nx, byrow = TRUE)[ny:1,]
            r <- raster::raster(mat, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=crs, template = NULL)
            return(r)
        },
        where = .rastermethods
    )

    methods::setMethod(raster::raster, methods::signature(x = "Dsurface"),
        ## define internal function for coercing Dsurface to raster
        function(x,  covariate, values = 1, crs = NA) {
            if (missing(covariate)) {
                covariate <- 'tempcov'
                covariates(x) <- data.frame(tempcov = rep(values, length.out=nrow(x)))
            }
            mask <- rectangularMask(x)
            bbox <- attr(x, 'boundingbox')
            xmx <- max(bbox$x)
            xmn <- min(bbox$x)
            ymx <- max(bbox$y)
            ymn <- min(bbox$y)
            nx <- round((xmx - xmn) / spacing(mask))
            ny <- round((ymx - ymn) / spacing(mask))
            tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
            x <- matrix (tmp, nrow = ny, ncol = nx, byrow = TRUE)[ny:1,]
            raster::raster(x, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=crs)
        },
        where = .rastermethods
    )
}

## checks
## getMethod('raster','mask')
## getMethod(raster, 'matrix')
## plot(raster(predictDsurface(fit.Dforest), 'D.0'), useRaster = F)

## 2022-01-24 terra SpatRaster
registerRastMethods <- function () {
    methods::setMethod(terra::rast, methods::signature(x = "mask"),
        ## define internal function for coercing mask to SpatRaster
        function(x,  covariate, values = 1, crs = "") {
            if (missing(covariate)) {
                covariate <- 'tempcov'
                covariates(x) <- data.frame(tempcov = rep(values, length.out=nrow(x)))
            }
            mask <- rectangularMask(x)
            tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
            df <- cbind(mask, tmp)
            terra::rast(df, type = 'xyz', crs = crs)
        },
        where = .rastermethods
    )

    methods::setMethod(terra::rast, methods::signature(x = "Dsurface"),
        ## define internal function for coercing Dsurface to SpatRaster
        function(x,  covariate, values = 1, crs = "") {
            if (missing(covariate)) {
                covariate <- 'tempcov'
                covariates(x) <- data.frame(tempcov = rep(values, length.out=nrow(x)))
            }
            mask <- rectangularMask(x)
            tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
            df <- cbind(mask, tmp)
            terra::rast(df, type = 'xyz', crs = crs)
        },
        where = .rastermethods
    )
}

## register now if the package is already loaded, otherwise when it is loaded
registerWhenLoaded <- function (pkg, fn) {
    if (isNamespaceLoaded(pkg)) {
        fn()
    }
    else {
        setHook(packageEvent(pkg, "onLoad"), function (...) fn())
    }
}
