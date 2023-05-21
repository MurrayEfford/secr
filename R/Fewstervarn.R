############################################################################################
## package 'secr'
## fewstervarn.R
## adapted in part from RMF
## simplified from derivedSystematic 2018-12-24-27
## 2019-12-22 lacework option
## 2019-12-23 single-detector cluster option
## 2019-12-29 ncores
############################################################################################

Fewstervarn <- function (nj, xy, design, esa, detectfn, detectpar, nocc,
                         basenx, df, maskspacing, keep, ncores, extrapolate = TRUE)
    {

    # design is a list with components --
    #   cluster -- traps object for one cluster
    #   region  -- sfc
    #   spacing -- distance between clusters)
    ## or --
    #   lacework -- logical TRUE
    #   region  -- sfc
    #   spacing -- 2-vector (distance between lines, between detectors along lines)
    # and possibly --
    #   exclude 
    #   edgemethod
    #   exclmethod
    design$region <- boundarytoSF(design$region)
    minxy <- st_bbox(design$region)[1:2]
    lacework <- !is.null(design$lacework) && design$lacework
    setNumThreads(ncores)  ## for multithreading in pdot
    #########################################################################
    ## Discrete sample space (bases for systematic samples)

    rx <- design$spacing[1]          ## spacing between cluster centres or lacework lines
    drx <- rx / basenx               ## spacing between base points
    bx <- seq(drx/2, rx, drx)
    b <- expand.grid(x = bx, y = bx) ## assume square grid
    B <- basenx^2                    ## number of base points
    
    ## 2019-12-26
    b <- b-rx

    #########################################################################
    ## Boxlet centres clipped to region
    
    boxlets <- make.mask(
        type = 'polygon', 
        poly = design$region, 
        spacing = maskspacing)
    ## optionally clip out 'lakes'
    if (!is.null(design$exclude)) {
        boxlets <- subset(boxlets, !pointsInPolygon(boxlets, design$exclude))
    }

    #########################################################################
    ## Model trend, using esa(j) to represent 'size' of each sample
    
    njxy <- data.frame (nj = nj, xy)  ## dataframe with one row per cluster
    if (is.null(df)) {
        lgam <- gam (nj ~ s(x, y), data = njxy,
                     family = poisson(link = log), offset = log(esa))
    }
    else {
        lgam <- gam (nj ~ s(x, y, fx = TRUE, k = df+1), data = njxy,
                     family = poisson(link = log), offset = log(esa))
    }
    
    #########################################################################
    ## Sum over boxlets for each basis point (starting position b)
    ## bxy is 2-vector of x- and y- displacement corresp to b
    ## pdot (X, trps, ...) gives p.(x) = 1 - (1-p(x))^nocc for animal at X
    ## where p(x) is aggregate over all detectors
    ## msk is subset of region near detectors

    sumoverboxlets <- function (bxy) {
        ## place grid for base 'bxy' and clip to region
        design$origin <- bxy + minxy 
        design$lacework <- NULL ## not an argument of make.lacework
        makefn <- if (lacework) 'make.lacework' else 'make.systematic'
        trps.b <- do.call(makefn, design)
        ## proceed cluster by cluster for speed, and to reduce memory demand
        onecluster <- function (msk, trps) {
            pxy <- covariates(msk)$pxy
            a <- pdot (msk, trps, detectfn, detectpar, nocc, ncores = ncores)
            c(sum(pxy * a), sum(a) * attr(msk, 'area'))
        }
        splittrps <- split(trps.b, clusterID(trps.b))
        splitbox  <- split(boxlets, clusters = splittrps)
        bycluster <- mapply(onecluster, splitbox, splittrps, SIMPLIFY = TRUE)
        ## return 2-vector of sum(pxy.a), sum(a) over clusters
        apply(bycluster, 1, sum)   
    }

    #########################################################################
    ## predicted p, saved as covariate of boxlet mask

    pxy <- predict.gam(lgam, newdata = boxlets, type = "response")
    covariates(boxlets) <- data.frame(pxy = pxy)
    if (!extrapolate) {
        inner <- xy[chull(xy),]
        inner <- pointsInPolygon(boxlets, inner)
        innerboxlets <- subset(boxlets, inner)
        outerboxlets <- subset(boxlets, !inner)
        nearest <- nearesttrap(outerboxlets, innerboxlets)
        covariates(boxlets)$pxy[!inner] <- covariates(boxlets)$pxy[nearest]
    }
    covariates(boxlets)$pxy <- covariates(boxlets)$pxy / 
        sum(covariates(boxlets)$pxy, na.rm=TRUE)
    
    #########################################################################
    ## Qb including pdot weight
    detectfn <- valid.detectfn(detectfn)  ## converts character code to numeric 
    tmp <- apply(b, 1, sumoverboxlets)    ## 2 x b matrix
    Qb <- tmp[1,]
    Ab <- tmp[2,]

    #########################################################################
    ## Finally...

    A <- sum(esa)
    regionarea <- masksize(boxlets)       ## hectares or km
    N <- sum(nj) / A  * regionarea        ## H-T estimate
    out <- 1/B * sum ((N * Qb * (1-Qb) + N^2 * Qb^2) / Ab^2) -
        sum((N * Qb / Ab) / B)^2
    out <- out * A^2                      ## scale var(n/A) to var(n)
    
    #########################################################################
    ## And optionally save lots of intermediate values
    
    if (keep) {
        attr(out, 'xy') <- xy
        attr(out, 'N') <- N
        attr(out, 'design') <- design
        attr(out, 'b') <- b
        attr(out, 'boxlets') <- boxlets  # pxy is covariate of mask
        attr(out, 'Qb') <- Qb
        attr(out, 'Ab') <- Ab
    }
    out  # var (n)
}
############################################################################

derivedSystematic <- function ( object, xy, design = list(),
    basenx = 10, df = 9, extrapolate = TRUE, alpha = 0.05, loginterval = TRUE, 
    independent.esa = FALSE, keep = FALSE, ncores = NULL) {
    ## 'design' should comprise, e.g.,
    ## list(cluster = hollowsquare, spacing = 600, region = possumarea)
    ## or
    ## list(lacework = TRUE, spacing = c(500,50), region = possumarea)
    
    warning("derivedSystematic is experimental in secr 4.6", call. = FALSE)
    
    if (!inherits (object, 'secr') | !is.null(object$groups))
        stop ("requires fitted secr model without groups")
    else if (ms(object)) {
        ## multisession, cluster = session
        if (('session' %in% object$vars) | ('Session' %in% object$vars) |
            (!is.null(object$sessioncov))) {
            stop ("derivedSystematic assumes detection model constant over sessions")
        }
        der <-  derived(object, se.esa = TRUE, se.D = FALSE, bycluster = FALSE)
        nj <- sapply(object$capthist, nrow)    ## animals per session
        maskspacing <- spacing(object$mask[[1]])
        detectpar <- detectpar(object)[[1]]
        nocc <- dim(object$capthist[[1]])[2]
        if (is.null(design$region)) {
            design$region <- attr(object$mask[[1]], 'polygon')
        }
    }
    else {
        ## single session, native clusters
        ## except that lacework treats each detector as a cluster
        if (is.null(clusterID(object$capthist))) {
            warning("clusters not defined; individual detectors as clusters")
            clusterID(traps(object$capthist)) <- 1:nrow(traps(object$capthist))
        }
        der <-  derived(object, se.esa = TRUE, se.D = FALSE, bycluster = TRUE)
        nj <- cluster.counts(object$capthist)    ## animals per cluster
        maskspacing <- spacing(object$mask)
        detectpar <- detectpar(object)
        nocc <- dim(object$capthist)[2]
        if (is.null(design$region)) {
            design$region <- attr(object$mask, 'polygon')
        }
    }
    der <- do.call(rbind, der)
    esa <- der[grepl('esa', rownames(der)),1]
    se.esa <- der[grepl('esa', rownames(der)),2]
    J <- length(nj)                              ## number of clusters
    
    if (is.null(design$region))
        stop ("region not found")
    
    ## convert matrix to sfc if not already
    design$region <- boundarytoSF(design$region)      ## see utility.R
    design$exclude <- boundarytoSF(design$exclude)    ## may be NULL

    #########
    ## var(n)
    varn <- Fewstervarn (nj, xy, design, esa, object$detectfn, detectpar, nocc,
                         basenx, df, maskspacing, keep, ncores)
    ##########
    ## var (A)
    A <- sum(esa)
    if (independent.esa)
        varA <- sum(se.esa^2)
    else
        varA <- J^2 * sum(esa/A * se.esa^2)
    
    ## unitary A - experimental
    lacework <- !is.null(design$lacework) && design$lacework
    if (lacework) {
        der <- derived(object, se.esa = TRUE, se.D = FALSE, bycluster = FALSE)
        A <- der['esa','estimate']
        varA <- der['esa','SE.estimate']^2
    }
    
    ##################
    ## assemble output
    n <- sum (nj)                          ## total animals
    if (lacework)
        D <- der['D','estimate']
    else 
        D <- n / A 
    
    varD <- D^2 * (varn/n^2 + varA/A^2)
    
    temp <- data.frame(row.names = c('esa','D'), 
                       estimate = c(A,D), 
                       SE.estimate = c(varA,varD)^0.5)
    temp <- add.cl(temp, alpha, loginterval)
    temp$CVn <- c(NA, varn^0.5/n)
    temp$CVa <- c(NA, sqrt(varA)/A)
    temp$CVD <- c(NA, varD^0.5/D)
    attr(temp, 'nj') <- nj
    attr(temp, 'esa') <- esa
    attr(temp, 'se.esa') <- se.esa
    
    ######################
    ## optional attributes
    if (keep) {
        attr(temp, 'xy') <- xy
        attr(temp, 'design') <- attr(varn, 'design')
        attr(temp, 'N') <- attr(varn, 'N')
        attr(temp, 'b') <- attr(varn, 'b')
        attr(temp, 'boxlets') <- attr(varn, 'boxlets')  # pxy is covariate of mask
        attr(temp, 'Qb') <- attr(varn, 'Qb')
        attr(temp, 'Ab') <- attr(varn, 'Ab')
    }
        
    temp
}

plotSystematic <- function (out, dec = 0, legend = TRUE, textcex = 0.6) {
    nj <- attr(out, 'nj')
    esa <- attr(out, 'esa')
    b <- attr(out, 'b')
    xy <- attr(out, 'xy')
    N <- attr(out, 'N')
    design <- attr(out, 'design')
    boxlets <- attr(out, 'boxlets')
    Qb <- attr(out, 'Qb')
    
    par(mfrow = c(2,3), mar = c(2,1,2,1), cex = 0.7)
    
    ##################################################
    ## nj
    plot(as(design$region, "Spatial"))  # avoid sf::plot
    text (xy[,1], xy[,2], as.character(nj), cex = textcex)
    mtext(side = 3, line = -0.5, 'nj')
    
    ##################################################
    ## nj/esa
    plot(as(design$region, "Spatial"))   # avoid sf::plot
    text (xy[,1], xy[,2], as.character(round(nj/esa,1)), cex = textcex)
    mtext(side = 3, line = -0.5, 'nj/esa')
    
    ##################################################
    ## base points b Not reliable for lacework
    opar <- par(mar=c(5,5,1,1))
    plot(as(design$region, "Spatial"))   # avoid sf::plot
    bbox <- st_bbox(design$region)
    points(b[,1] + bbox[1], 
           b[,2] + bbox[2], xpd = TRUE)
    mtext(side = 3, line = -0.5, 'base points b')
    par(opar)
    
    ##################################################
    ## grid for b[1]
    design$origin <- b[1,] + bbox[1:2]
    lacework <- !is.null(design$lacework) && design$lacework
    design$lacework <- NULL
    if (lacework) {
        trps <- do.call(make.lacework, design)
    }
    else
        trps <- do.call(make.systematic, design)
    plot(as(design$region, "Spatial"))   # avoid sf::plot
    plot(trps, detpar = list(cex=0.7), add = TRUE)
    mtext(side = 3, line = -0.5, 'grid for b[1]')
    
    ##################################################
    ## pxy
    covariates(boxlets)$pxy4 <- covariates(boxlets)$pxy * 1e4
    plot(boxlets)
    plot(boxlets, cov = 'pxy4', breaks = 5, add = TRUE, legend = legend)
    mtext(side = 3, line = -0.5, 'pxy * 1E4')
    plot(as(design$region, "Spatial"), add = TRUE)   # avoid sf::plot

    ##################################################
    ## E(n|b)
    opar <- par(pty = 's', mar = c(4,4,2,1))
    plot(b, type = 'n', xlab = '', ylab = '', axes = FALSE, 
         xlim = c(-design$spacing[1],0), ylim=c(-design$spacing[1],0))
    axis(1, at = unique(b$x))
    axis(2, at = unique(b$y))
    text (b[,1], b[,2], round(N*Qb, dec+1), cex = 0.7)
    mtext(side = 3, line = 1, 'E(n|b)')
    par(opar)
    ##################################################
    
}

# plotSystematic(tmp)

