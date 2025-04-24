## Simulated categorical landscape
## Modified Random Cluster method of Saura and Martinez-Millan 2000
## Landscape Ecology 15: 661-678

## Murray Efford 2012-04-09,10,11,12;
## 2013-07-09 replaced call to obsolete raster function 'adjacency'
## 2014-08-25 'raster' now in Imports, so do not 'require'
## 2017-07-26 seed argument suggested by Erin Peterson

## Requires 'raster' and 'igraph0' packages:

## Robert J. Hijmans & Jacob van Etten (2011). raster: Geographic
##  analysis and modeling with raster data. R package version 1.9-33.
##  https://CRAN.R-project.org/package=raster

## Csardi G, Nepusz T: The igraph software package for complex network
##  research, InterJournal, Complex Systems 1695. 2006.
##  https://igraph.org

## Restrictions
##   adjacency for step 'B' is based on a 1-step rook move (direction = 4)
##   (this yields an isotropic result)

## 2023-08-23 randomDensity

randomHabitat <- function (mask, p = 0.5, A = 0.5, directions = 4, minpatch = 1,
                     drop = TRUE, covname = 'habitat', plt = FALSE, seed = NULL) {

    #############################
    ## set random seed 2017-07-26
    ## copied from simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    #############################
    
    if (ms(mask)) {
        ## allow for possibility that 'mask' is a list of masks
        temp <- lapply(mask, randomHabitat, p = p, A = A, directions = directions,
                       minpatch = minpatch, drop = drop, covname = covname, plt = plt,
                       seed = NULL)
        class(temp) <- c('mask', 'list')
        attr(temp,'seed') <- RNGstate   ## save random seed
        temp
    }
    else {

        if (abs(A-1) < 1e-6) {
            ## 2023-08-23 allow A=1 for 100% habitat
            covariates(mask)[[covname]] <- rep(1, nrow(mask))    
            return(mask)
        }
        
        ## extract limits etc. from input mask
        spacing <- attr(mask,'area')^0.5 * 100
        bb <- attr(mask, 'boundingbox')
        maxx <- max(bb$x)
        minx <- min(bb$x)
        maxy <- max(bb$y)
        miny <- min(bb$y)
        nx <- round( (maxx - minx) / spacing )
        ny <- round( (maxy - miny) / spacing )
        n <- nx * ny

        ## create rasterLayer
        layer <- raster(nrows = ny, ncols = nx, xmn = minx, xmx = maxx,
                                ymn = miny, ymx = maxy)

        ## A. Percolation map generation
        raster::values(layer) <- rep(0, n)
        raster::values(layer)[sample.int(n, round(n*p))] <- 1   ## as close as possible
        if (plt) raster::plot(layer, useRaster = FALSE)

        ## B. Cluster identification (single-linkage clustering of adjoining pixels)
        if (!requireNamespace("igraph", quietly = TRUE)) {
            stop("you need to install the igraph package to use randomHabitat")
        }
        clumped <- raster::clump(layer, directions = directions, gaps = FALSE)

        ## C. Cluster type assignment
        ncluster <- max(raster::values(clumped), na.rm = TRUE)
        types <- factor(c(0,1))          # nonhabitat, habitat
        numTypes <- as.numeric(types)    # 1,2
        clustertype <- sample(numTypes, ncluster, replace = TRUE, prob = c(1-A,A))
        raster::values(clumped) <- clustertype[raster::values(clumped)]
        if (plt) raster::plot(clumped, useRaster = FALSE)

        ## D. Filling in image
        cellsUnassigned <- (1:n)[is.na(raster::values(clumped))]
        cellsAssigned <- (1:n)[!is.na(raster::values(clumped))]
        tempadj <- raster::adjacent (clumped, cells = cellsUnassigned, target = cellsAssigned,
                             directions = 8)
        tempadj <- split(tempadj[,2], tempadj[,1])
        fillinType <- function (adjcells) {
            type <- raster::values(clumped)[adjcells]
            type <- factor(type)
            freq <- tabulate(type)
            result <- as.numeric(levels(type)[freq == max(freq)])
            if (length(result)>1)
                result <- sample (result, 1)
            result
        }
        # cells with typed neighbours
        filled <- sapply(tempadj, fillinType)
        filledCells <- as.numeric(names(filled))
        raster::values(clumped)[filledCells] <- filled
        # cells with no typed neighbours
        notfilledCells <- cellsUnassigned[!(cellsUnassigned %in% filledCells)]
        randomType <- sample(numTypes, length(notfilledCells), replace = TRUE, prob = c(1-A,A))
        raster::values(clumped)[notfilledCells] <- randomType
        if (plt) {
            # flip to match secr plotting (y=0 at bottom) 2023-08-23
            raster::plot(raster::flip(clumped), useRaster = FALSE)
        }
        raster::values(clumped) <- raster::values(clumped) - 1

        ## optionally filter small patches
        if (minpatch > 1) {
            reclumped <- raster::clump(clumped, directions = directions, gaps = FALSE)
            nperclump <- table(raster::values(reclumped))
            raster::values(clumped)[nperclump[as.character(raster::values(reclumped))] < minpatch] <- 0
            temp <- clumped
            raster::values(temp) <- 1-raster::values(temp)   ## swap 0,1
            reclumped <- raster::clump(temp   , directions = directions, gaps = FALSE)
            nperclump <- table(raster::values(reclumped))
            raster::values(clumped)[nperclump[as.character(raster::values(reclumped))] < minpatch] <- 1
        }

        # pad if necessary (sets attribute OK with padding cells FALSE)
        rectmask <- rectangularMask(mask)
        # restrict to 'habitat'; discard padding
        tempmask <- subset(rectmask, attr(rectmask, 'OK') & (raster::values(clumped)>0))
        # optionally return mask with only 'habitat' points
        if (drop) {
            mask <- tempmask
            attr(mask, 'type') <- paste('MRC p=',p, ' A=',A, sep='')
        }
        else {
            hab <-  as.numeric(pointsInPolygon (mask, tempmask))
            if (is.null(covariates(mask))) {
                covariates(mask) <- data.frame(hab)
                names(covariates(mask)) <- covname
            }
            else
                covariates(mask)[,covname] <- hab
        }
        attr(mask, 'seed') <- RNGstate   ## save random seed
        
        mask
    }
}
################################################################################

randomDensity <- function (mask, parm) 
{
    defaultparm <- list(D = NULL, p = 0.5, A = 0.5, directions = 4, 
                        minpatch = 1, plt = FALSE, seed = NULL, 
                        rescale = TRUE)
    if (is.null(parm$D)) stop ("randomDensity requires D to be specified")
    parm     <- replace(defaultparm, names(parm), parm)
    userargs <- c("p", "A", "directions", "minpatch", "plt", "seed")
    tempmask <- do.call(randomHabitat, 
                        c(parm[userargs], mask = list(mask), drop = FALSE))
    habitat  <- covariates(tempmask)[["habitat"]]
    if (parm$rescale)
        # adjust for requested habitat proportion
        habitat * parm$D / parm$A
    else
        habitat * parm$D
}
################################################################################


