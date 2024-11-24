##############################################################################
## package 'secr'
## addCovariates.R
## 2011-11-01
## 2013-01-19 handles missing values in character fields
## 2014-08-05 strict argument
## 2014-08-05 relax requirement for object to be traps or mask:
## 2014-08-05 may now be any numeric vector that can be formed into a 2-column matrix
## 2017-03 new argument replace; use readOGR for shapefiles
## 2021-09-16 allow raster
## 2021-12-07 SpatRaster (terra)
## 2022-02-13 Revamped for sf, added tests
## 2022-08-28 extended to popn objects
## 2024-11-13 fixed st_join() call to avoid spurious warning "attribute variables 
#             are assumed to be spatially constant throughout all geometries"
###############################################################################

addCovariates <- function (object, spatialdata, columns = NULL, strict = FALSE, replace = FALSE) {

    if (!(inherits(object, c('mask', 'traps', 'popn')))) 
        object <- matrix(unlist(object), ncol = 2)
    if (!ms(object) & ms(spatialdata))
        stop ("mismatch of single session object, multisession spatialdata")

    #---------------------------------------------------------------------------
    # multisession
    if (ms(object)) {
        ## allow multiple sessions, and session-specific data sources
        nsession <- length(object)
        out <- object
        for (session in 1:nsession) {
            if (ms(spatialdata)) {
                out[[session]] <- addCovariates(out[[session]], spatialdata[[session]], 
                    columns, strict, replace)
            }
            else {
                out[[session]] <- addCovariates(out[[session]], spatialdata, 
                    columns, strict, replace)
            }
        }
        out
    }
    #---------------------------------------------------------------------------
    # single session
    else {
        if (inherits(spatialdata, c('SpatRaster','Rasterlayer', 'SpatialGridDataFrame'))) {
            if (!requireNamespace('terra', quietly = TRUE)) {
                stop ("package 'terra' >= 1.5-12 is required to add covariates from a raster data source")
            }
        }
        # transform spatial data to sf, SpatRaster, traps or mask
        if (is.character(spatialdata)) {
            polyfilename <- spatialdata  
            isshp <- function(filename) {
                nch <- nchar(filename)
                tolower(substring(filename, nch-3,nch)) == ".shp" 
            }
            if (!isshp(polyfilename)) {
                polyfilename <- paste0(polyfilename, ".shp")
            }
            spatialdata <- st_read(polyfilename, quiet = TRUE)   # read sf
        }
        else if (inherits(spatialdata, "SpatialPolygonsDataFrame")) {
            spatialdata <- st_as_sf(spatialdata)
        }
        else if (inherits(spatialdata, "SpatialGridDataFrame")) {
            spatialdata <- terra::rast(raster(spatialdata))
        }
        else if (inherits(spatialdata, "RasterLayer")) {
            spatialdata <- terra::rast(spatialdata)
        }
        # process each allowed type
        if (inherits(spatialdata, "sf")) {
            # POLYGON or MULTIPOLYGON
            xy <- as.data.frame(object)
            xy <- st_as_sf(xy, coords=1:2, crs = st_crs(spatialdata))
            # removed 'largest' argument 2024-11-13
            df <- st_join(xy, spatialdata, join = st_within)
            df <- st_drop_geometry(df)
        }
        else if (inherits(spatialdata, "SpatRaster")) {
            df <- data.frame(raster = terra::extract(spatialdata, as.matrix(object)))
            if (!is.null(columns)) {
                names(df) <- columns
            }
        }
        else if (inherits(spatialdata, c("traps", "mask"))) {
            ## nearest point algorithm
            if (is.null(covariates(spatialdata)))
                stop ("spatialdata does not have covariates")
            index <- nearesttrap(object, spatialdata)
            df <- covariates(spatialdata)[index,, drop=FALSE]
            ## new argument 2014-08-05
            if (strict && inherits(spatialdata, "mask")) {
                incell <- function (xy, m, mask) {
                    sp2 <- spacing(mask) / 2
                    mxy <- mask[m,]
                    ((xy[,1] + sp2) >= mxy[,1]) &
                    ((xy[,1] - sp2) <= mxy[,1]) &
                    ((xy[,2] + sp2) >= mxy[,2]) &
                    ((xy[,2] - sp2) <= mxy[,2])

                }
                cellOK <- incell(object, index, spatialdata)
                df[!cellOK,] <- NA
                if (any(!cellOK))
                    warning ("some requested points lie outside mask")
            }
        }
        else {
            stop ("spatialdata type unrecognised or unsupported")
        }
        
        ## select requested columns
        if (!is.null(columns)) {
            df <- df[,columns, drop = FALSE]
        }

        ## check new covariates OK
        fn <- function(x) {
            if (is.numeric(x))
                !any(is.na(x))
            else
                !any((x == "") | is.na(x))
        }
        OK <- all(apply(df, 2, fn))
        if (!OK) {
            warning ("missing values among new covariates")
        }
        
        ## insert new covariates and return object
        rownames(df) <- 1:nrow(df)
        if (is.null(covariates(object)))
            covariates(object) <- df
        else {
            if (replace) {
                repeated <- names(covariates(object)) %in% names(df)
                covariates(object) <- covariates(object)[,!repeated]
            }
            
            covariates(object) <- cbind(covariates(object), df)
        }
        object
    }
}

