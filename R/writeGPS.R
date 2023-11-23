###############################################################################
## package 'secr'
## writeGPS.R
## last changed 2022-02-01 (sf)
## Write detector locations to GPS format
##
###############################################################################

## writeGPS basics
## o = 'garmin' is only for Garmin usb/serial protocol
## use o = 'gdb' for MapSource

writeGPS <- function (xy, o = "garmin", F = "usb:", proj = '+proj=nzmg')  {

    if (proj=="") {
        latlon <- xy
    }
    else {
        ## express as latitude/longitude
        xy <- sf::st_as_sf(as.data.frame(xy), coords = 1:2, crs = proj)
        latlon <- sf::st_transform(xy, "epsg:4326")
        latlon <- sf::st_coordinates(latlon)
    }
    latlon <- data.frame(latlon)[,c(2,1)] ## need lat first
    latlon$label <- rownames(xy)
    tempf <- tempfile('waypts')
    oldopt <- options(digits=12)  ## ensure plenty of digits
    write.table(latlon, quote = FALSE, file = tempf,
        col.names = FALSE, row.names = FALSE, sep = ',')

    ## adapted from maptools::readGPS
    GB <- Sys.which("gpsbabel")
    if (nchar(GB) == 0 || !file.exists(GB))
        stop ("gpsbabel not found")
    cmd <- paste(GB, "-D9 -w -i csv -f ", tempf, " -o ", o, " -F ", F)
    if (.Platform$OS.type == "windows")
        gpsdata <- system(cmd, intern = TRUE, invisible = TRUE)
    else
        gpsdata <- system(cmd, intern = TRUE)

    file.remove(tempf)
    options(oldopt)
    invisible()

}
