#######################################################################################
## package 'secr'
## gridcells.R
## 2022-05-30
#######################################################################################

gridCells <- function (x, cellsize = spacing(x), crs = NA) {
    dx <- cellsize / 2
    dy <- cellsize / 2
    makepixel <- function (i) {
        xy <- unlist(x[i,])
        mat <- cbind(
            x = xy[1] + dx*c(-1,-1,1,1,-1),
            y = xy[2] + dy*c(-1,1,1,-1,-1))
        sf::st_linestring(mat)
    }
    pixels <- lapply(1:nrow(x),makepixel)
    pixellines <- do.call(sf::st_sfc, pixels)
    pixellines <- sf::st_cast(pixellines, 'MULTILINESTRING')
    cells <- sf::st_cast(pixellines, 'MULTIPOLYGON')
    sf::st_crs(cells) <- crs
    cells
}
