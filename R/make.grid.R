###############################################################################
## package 'secr'
## make.grid.R
## make.grid, make.poly, make.transect, make.circle
## 2016-10-16 secr 3.0
## 2017-05-23 markocc argument
## 2019-12-29 leadingzero argument
## 2021-03-26 spacing attribute initially inferred from min(spacex, spacey)
###############################################################################

make.grid <- function (nx = 6, ny = 6, spacex = 20, spacey = spacex, spacing=NULL, 
    detector = 'multi', originxy = c(0,0), hollow = FALSE, ID = 'alphay', 
    leadingzero = TRUE, markocc = NULL)

{
    if (!all( detector %in% .localstuff$validdetectors ))
        stop ("invalid detector type")
    if (!is.null(spacing)) {
        spacex <- spacing
        spacey <- spacing
    }
    originxy <- as.numeric(unlist(originxy))  ## 2020-01-02 for robustness
    if ((nx<3) | (ny<3))
        hollow <- FALSE
    grid <- expand.grid (x=(0:(nx-1))*spacex + originxy[1], y=(0:(ny-1))*spacey + originxy[2])

    allowedID <- c('numx','numy', 'numxb', 'numyb', 'alphax','alphay', 'xy')
    if (! ID %in% allowedID)
        stop ("ID should be one of ",
            paste(sapply(allowedID, dQuote),collapse=','))

    if (ID %in% c('alphax','alphay')) {
        n <- max(nx,ny)/25 + 1
        ll <- sapply(1:n, function(x) apply(matrix(rep(LETTERS,rep(x,26)), nrow = x), 2,
            paste, collapse=''))
    }

    if (ID == 'numy') temp <- 1:nrow(grid)
    if (ID == 'numx') temp <- t(matrix(1:nrow(grid), ncol = nx))

    if (ID == 'numyb') {
        temp <- matrix(1:(nx*ny), ncol = ny)
        if (ny>1)  # 2012-07-13
        for (i in seq(2,ny,2)) temp[,i] <- rev(temp[,i])
    }

    if (ID == 'numxb') {  ## YES
        temp <- t(matrix(1:(nx*ny), ncol = nx))
        if (nx>1)  # 2012-07-13
            for (i in seq(2,nx,2)) temp[i,] <- rev(temp[i,])
    }

    if (ID == 'alphax') {
        colA <- ll[1:nx]
        rowA <- 1:ny
        row.names(grid) <- apply(expand.grid(colA,rowA), 1, paste, sep='', collapse='')
    }
    if (ID == 'alphay') {
        colA <- 1:nx
        rowA <- ll[1:ny]
        row.names(grid) <- apply(expand.grid(colA, rowA), 1,
            function(x) paste(rev(x), sep='', collapse=''))
    }

    ## added 2010 03 24
    if (ID == 'xy') {
        ## see utility.R for function secr_leadingzero 2012-09-04
        colA <- secr_leadingzero(1:nx)
        rowA <- secr_leadingzero(1:ny)
        row.names(grid) <- apply(expand.grid(colA, rowA), 1,
            function(x) paste(x, sep='', collapse=''))
    }

    if (hollow) {
        grid <- grid[grid$x==originxy[1] |
                             grid$x==(originxy[1] + spacex*(nx-1)) |
                             grid$y==originxy[2] |
                             grid$y==(originxy[2] + spacey*(ny-1)),]

        ## number clockwise from bottom left
        grid <- grid[order(grid$x, grid$y),]
        temp <- c(1:ny,
            t(matrix(c(
                rev ((2*ny+nx-1):(2*ny+2*nx-4)),
                (ny+1):(ny+nx-2)
              ), ncol = 2)),
            rev((ny+nx-1):(2*ny+nx-2)))
    }

    if (hollow | (ID %in% c('numx','numy','numxb','numyb'))) {
        row.names(grid) <- temp
        grid <- grid[order(temp),]
    }

    if (leadingzero && (ID %in% c('numy','numx','numyb','numxb') )) {
        row.names(grid) <- secr_leadingzero(as.numeric(row.names(grid)))
    }
    
    attr(grid, 'detector')    <- detector
    attr(grid, 'class')       <- c('traps', 'data.frame')
    attr(grid, 'spacex')      <- spacex
    attr(grid, 'spacey')      <- spacey
    attr(grid, 'spacing')     <- min(spacex, spacey)  # 2021-03-26
    attr(grid, 'spacing')     <- spacing(grid)  ## reset if NULL
    attr(grid, 'usage')       <- NULL
    attr(grid, 'cluster')     <- NULL
    attr(grid, 'clustertrap') <- NULL
    attr(grid, 'covariates')  <- NULL
    
    ##2016-10-05 to fix example in trapsfn.Rd
    if (any(detector %in% .localstuff$polydetectors))
        polyID(grid) <- rep(1, nrow(grid))

    if (!is.null(markocc))
        markocc(grid) <- markocc
    grid
}
###############################################################################

make.poly <- function (polylist=NULL, x=c(-50,-50,50,50), y=c(-50,50,50,-50),
                       exclusive = FALSE, verify = TRUE)
# polygon detectors
{
    makepart <- function (vert) {
        ## as.matrix to avoid as.data.frame.traps
        vert <- data.frame(as.matrix(vert))
        if ((ncol(vert)==2) & !all(c('x','y') %in% names(vert)))
            names(vert) <- c('x','y')
        if (any(tail(vert,1) != vert[1,]))   ## close polygon
            vert <- rbind(vert, vert[1,,drop = FALSE])
        vert
    }

    if (is.null(polylist)) polylist <- list(data.frame(x=x,y=y))
    if (is.null(names(polylist))) names(polylist) <- 1:length(polylist)
    grid <- lapply (polylist, makepart)
    polyn <- names(polylist)
    nrg <- unlist(sapply(grid, nrow))
    grid <- data.frame(abind(grid,along=1), row.names=NULL)
    class(grid)     <- c('traps', 'data.frame')
    if (exclusive)
        detector(grid)  <- 'polygonX'
    else
        detector(grid)  <- 'polygon'
    ## 2011-06-23
    ## polyID(grid) <- factor(rep(polyn, nrg))
    polyID(grid) <- factor(rep(polyn, nrg), levels = polyn)
    if (verify) verify(grid)
    grid
}
###############################################################################

make.telemetry <- function (xy = c(0,0)) {
    temp <- data.frame(x = xy[1], y = xy[2])
    class(temp)     <- c('traps', 'data.frame')
    detector(temp) <- 'telemetry'
    attr(temp, 'telemetrytype') <- "independent"
    attr(temp, 'spacex')      <- NULL
    attr(temp, 'spacey')      <- NULL
    attr(temp, 'spacing')     <- NULL
    attr(temp, 'usage')       <- NULL
    attr(temp, 'cluster')     <- NULL
    attr(temp, 'clustertrap') <- NULL
    attr(temp, 'covariates')  <- NULL
    temp
}

make.transect <- function (transectlist=NULL, x=c(-50,-50,50,50), y=c(-50,50,50,-50),
                           exclusive = FALSE)
# transect detectors
{
    makepart <- function (vert) {
        vert <- data.frame(vert)
        if ((ncol(vert)==2) & !all(c('x','y') %in% names(vert)))
            names(vert) <- c('x','y')
        vert
    }
    if (is.null(transectlist)) transectlist <- list(data.frame(x=x,y=y))
    if (is.null(names(transectlist))) names(transectlist) <- 1:length(transectlist)
    grid <- lapply (transectlist, makepart)
    transectn <- names(transectlist)
    nrg <- unlist(sapply(grid, nrow))
    grid <- data.frame(abind(grid,along=1), row.names=NULL)
    class(grid)     <- c('traps', 'data.frame')
    if (exclusive)
        detector(grid)  <- 'transectX'
    else
        detector(grid)  <- 'transect'
    ## 2011-06-23
    ## polyID(grid) <- factor(rep(transectn, nrg))
    polyID(grid) <- factor(rep(transectn, nrg), levels = transectn)
    grid
}
###############################################################################

make.circle <- function (n = 20, radius = 100, spacing = NULL,
    detector = 'multi', originxy=c(0,0), IDclockwise = TRUE, 
    leadingzero = TRUE, markocc = NULL)
{
    if (!all( detector %in% .localstuff$validdetectors ))
        stop ("invalid detector type")
    if (is.null(radius) & is.null(spacing))
        stop ("specify 'radius' or 'spacing'")
    theta <- seq (0, 2 * pi * (n-1)/n, 2 * pi / n)
    if (!is.null(spacing)) radius <- spacing / 2 / sin(pi / n)  ## override
    if (is.null(spacing)) spacing <- radius * sin(pi / n) * 2

    object <- data.frame (x = radius * cos(theta) + originxy[1],
                          y = radius * sin(theta) + originxy[2])

    if (IDclockwise) sequence  <- c(1, n:2)
    else sequence <- 1:n
    if (leadingzero)
        row.names(object) <- secr_leadingzero(sequence)
    else 
        row.names(object) <- sequence
    
    object <- object[order(sequence),]
    attr(object, 'detector')    <- detector
    attr(object, 'class')       <- c('traps', 'data.frame')
    attr(object, 'spacex')      <- spacing
    attr(object, 'spacey')      <- spacing
    attr(object, 'spacing')     <- spacing
    attr(object, 'usage')       <- NULL
    attr(object, 'cluster')     <- NULL
    attr(object, 'clustertrap') <- NULL
    attr(object, 'covariates')  <- NULL
    if (!is.null(markocc))
        markocc(object) <- markocc
    object
}
###############################################################################


