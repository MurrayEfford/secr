############################################################################################
## package 'secr'
## dbar.R, ARL.R, MMDM.R, RPSV.R
## all in homerange.R 2014-09-01
## 2014-09-01 modified for userdist
## 2014-09-10 does NOT work when userdist fn requires mask covariates..
## 2016-10-06 secr 3.0
## 2016-10-28 NOT YET userdist may be session-specific
## 2017-02-06 updated for telemetry (completing job started last month)
## 2020-08-31 ORL, centroids
## 2021-05-18 sortorder
## 2025-06-16 t2r2
############################################################################################

getID <- function (det, capthist) {
    if (all(det %in% .localstuff$polydetectors)) {
        ID <- animalID(capthist, names = TRUE, sortorder = "ksn")
    }
    else {
        ID <- animalID(capthist, names = TRUE, sortorder = "snk")
    }
    factor(ID, levels = rownames(capthist))
}

dbar <- function (capthist, userdist = NULL, mask = NULL) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, dbar, userdist, mask)   ## recursive
    }
    else {
        dbarx    <- function (x) {
            x <- abs(unlist(x))
            ## assume that within animal, x is in order by occasion
            distmat[cbind(x[-length(x)], x[-1])]  ## vector
        }
        dbarxy    <- function (xy) {
            sqrt(diff(xy$x)^2 + diff(xy$y)^2)
        }
        if (nrow(capthist) < 1) return(NA)
        traps <- traps(capthist)
        det <- secr_expanddet(capthist)
        IDfactor <- getID(det, capthist)
        ## 2014-09-01
        ## NOT USING PARAMETERS noneuc ETC
        distmat <- secr_valid.userdist(userdist, det, traps, traps, mask )
        if (!all(det %in% .localstuff$individualdetectors))
            stop ("require individual detector type for dbar")
        
        if (all(det %in% 'telemetry')) {
            lxy <- telemetryxy(capthist)
            if (is.null(lxy))
                NA
            else {
                d <- try(lapply(lxy,dbarxy), silent = TRUE)
                if (inherits(d, 'try-error'))
                    d <- NA
                mean(unlist(d), na.rm=T)
            }
        }
        else if (all(det %in% .localstuff$polydetectors)) {
            if (is.null(xy(capthist)))
                NA
            else {
                lxy <- split (xy(capthist), IDfactor)
                d <- try(lapply(lxy,dbarxy), silent = TRUE)
                if (inherits(d, 'try-error'))
                    d <- NA
                mean(unlist(d), na.rm=T)
            }
        }
        else {
            ## order is essential 2016-10-07
            if (any(det %in% 'telemetry')) {
                capthist <- subset(capthist, 
                                   occasions = det != 'telemetry', 
                                   traps = 1:(nrow(traps(capthist)-1)))
                IDfactor <- factor(animalID(capthist, names = TRUE), 
                    levels = rownames(capthist))
            }
            w <- split(trap(capthist, names = FALSE), IDfactor)
            d <- try(unlist(lapply(w,dbarx)), silent = TRUE)
            if (inherits(d, 'try-error'))
                d <- NA
            mean(d, na.rm=T)
        }
    }
}
############################################################################################

moves <- function (capthist, userdist = NULL, mask = NULL, names = FALSE) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, moves, userdist, mask)   ## recursive
    }
    else {
        movex    <- function (x) {
            x <- abs(unlist(x))
            distmat[cbind(x[-length(x)], x[-1])]  ## vector
        }
        movexy    <- function (xy) {
            sqrt(diff(xy$x)^2 + diff(xy$y)^2)
        }
        traps <- traps(capthist)
        if (is.null(traps)) {
            ## return empty vector if nonspatial 2019-04-04
            nam <- rownames(capthist)
            if (is.null(nam)) nam <- 1:nrow(capthist)
            sapply(nam, function(x) numeric(0), simplify = FALSE)
        }
        else {
            det <- secr_expanddet(capthist)
            IDfactor <- getID(det, capthist)
            distmat <- secr_valid.userdist(userdist, det, traps, traps, mask)
            if (!all(det %in% .localstuff$individualdetectors))
                stop ("require individual detector type for moves")
            
            if (all(det %in% 'telemetry')) {
                lxy <- telemetryxy(capthist)
                if (is.null(lxy))
                    out <- NA
                else {
                    out <- lapply (lxy, movexy)
                }
            }
            else if (all(det %in% .localstuff$polydetectors)) {
                if (is.null(xy(capthist)))
                    out <- NA
                else {
                    lxy <- split (xy(capthist), IDfactor)
                    out <- lapply (lxy, movexy)
                }
            }
            else {
                ## order is essential 2016-10-08
                if (any(det %in% 'telemetry')) {
                    capthist <- subset(capthist, 
                                       occasions = det != 'telemetry', 
                                       traps = 1:(nrow(traps(capthist)-1)))
                    IDfactor <- factor(animalID(capthist, names = TRUE), 
                        levels = rownames(capthist))
                }
                ## 2020-08-27
                w <- split(trap(capthist, names = FALSE), IDfactor)
                out <- lapply(w,movex)
            }
            ## 2022-01-21
            # if (!names) names(out) <- 1:length(out)
            if (!names && length(out)>0) names(out) <- 1:length(out)
            out
        }
    }
}
############################################################################################

trapsPerAnimal <- function (capthist) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, trapsPerAnimal)   ## recursive
    }
    else {
        nki <- apply( apply(abs(capthist), c(1,3), sum)>0, 1, sum)
        if (length(nki)>0)
            out <- tabulate (nki, nbins = max(nki))
        else
            out <- 0
        names(out) <- 1:length(out)
        out
    }
}
############################################################################################


ARL <- function (capthist, min.recapt = 1, plt = FALSE, full = FALSE, userdist = NULL,
                 mask = NULL) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, ARL, plt = plt, full = full, userdist = userdist, mask)   ## recursive
    }
    else {
        MMDMx <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            if (sum(cx>0, na.rm=T) == 1) NA
            else {
              ## x <- traps$x[cx]
              ## y <- traps$y[cx]
              ## max(dist(cbind(x,y)))
              as.numeric(max(distmat[cx, cx]))
            }
        }
        MMDMxy <- function (xy) {
            max(dist(cbind(xy$x, xy$y)))
        }
        if (nrow(capthist) < 1) return(NA)
        traps <- traps(capthist)
        det <- secr_expanddet(capthist)
        IDfactor <- getID(det, capthist)
        if (!all(det %in% .localstuff$individualdetectors))
            stop ("require individual detector type for ARL")
        distmat <- secr_valid.userdist(userdist, det, traps, traps, mask )
        prox  <- length(dim(capthist)) > 2
        
        if (all(det %in% 'telemetry')) {
            lxy <- telemetryxy(capthist)
            if (is.null(lxy))
                stop("no telemetry coordinates")
            else {
                maxd <- unlist(lapply (lxy, MMDMxy))
                n <- unlist(lapply (lxy, nrow))
            }
        }
        else if (all(det %in% .localstuff$polydetectors)) {
            if (is.null(xy(capthist)))
                stop("no xy coordinates")
            else {
                lxy <- split (xy(capthist), IDfactor)
                maxd <- unlist(lapply (lxy, MMDMxy))
                n <- unlist(lapply (lxy, nrow))
            }
        }
        else {
            ## order is essential 2016-10-08
            if (any(det %in% 'telemetry')) {
                capthist <- subset(capthist, 
                                   occasions = det != 'telemetry', 
                                   traps = 1:(nrow(traps(capthist)-1)))
                IDfactor <- factor(animalID(capthist, names = TRUE), 
                    levels = rownames(capthist))
            }
            w <- split(trap(capthist, names = FALSE), IDfactor)
            maxd <- unlist(lapply(w, MMDMx))
            n <- unlist(lapply(w, length))
        }
        maxd <- maxd[n>min.recapt]
        n <- n[n>min.recapt]

        temp <- try(coef(nls (maxd ~ aa * (1 - exp(bb * (n-1))),
            start= list (aa = max(maxd)*1.2, bb = -0.4))))

        if (inherits(temp, "try-error")) {
            warning ("nls failure")
            aa <- NA
            bb <- NA
        }
        else {
            aa <- temp[1]
            bb <- temp[2]
            if (plt) {
                plot (jitter(n, amount=0.1), maxd,
                    xlim=c(0,max(c(n,ncol(capthist)))),
                    xlab='Number of captures', ylab='ARL')
                xv <- seq(2,max(n),0.01)
                lines(xv, aa * (1 - exp(bb * (xv-1))))
            }
        }
        attr(aa,'names') <- NULL
        attr(bb,'names') <- NULL

        if (!full) aa
        else list (ARL = aa, b = bb, data = data.frame(maxd = maxd, n=n))
    }
}
############################################################################################

MMDM <- function (capthist, min.recapt = 1, full = FALSE, userdist = NULL, mask = NULL) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, MMDM, full = full, userdist = userdist, mask = mask)   ## recursive
    }
    else {
        MMDMx <- function (cx) {
              cx <- abs(cx)  # no special trt for deads
              if (sum(cx>0, na.rm=T) == 1) NA
              else {
                ## x <- traps$x[cx]
                ## y <- traps$y[cx]
                ## max(dist(cbind(x,y)))
                as.numeric(max(distmat[cx, cx]))
              }
        }
        MMDMxy    <- function (xy) {
                max(dist(cbind(xy$x, xy$y)))
        }
        if (nrow(capthist) < 1) return(NA)
        traps <- traps(capthist)
        det <- secr_expanddet(capthist)
        IDfactor <- getID(det, capthist)
        distmat <- secr_valid.userdist(userdist, det, traps, traps, mask )
        if (!all(det %in% .localstuff$individualdetectors))
            stop ("require individual detector type for MMDM")
        
        if (all(det %in% 'telemetry')) {
            lxy <- telemetryxy(capthist)
            if (is.null(lxy))
                stop ("no telemetry coordinates")
            else {
                maxd <- unlist(lapply (lxy, MMDMxy))
                n <- unlist(lapply (lxy, nrow))
            }
        }
        else if (all(det %in% .localstuff$polydetectors)) {
            if (is.null(xy(capthist)))
                stop ("no xy coordinates")
            else {
                lxy <- split (xy(capthist), IDfactor)
                maxd <- unlist(lapply (lxy, MMDMxy))
                n <- unlist(lapply (lxy, nrow))
            }
        }
        else {
            if (any(det %in% 'telemetry')) {
                capthist <- subset(capthist, 
                                   occasions = det != 'telemetry', 
                                   traps = 1:(nrow(traps(capthist)-1)))
                IDfactor <- factor(animalID(capthist, names = TRUE), 
                    levels = rownames(capthist))
            }
            w <- split(trap(capthist, names = FALSE), IDfactor)
            maxd <- unlist(lapply( w, MMDMx))
            n <- unlist(lapply(w, length))
        }
        maxd <- maxd[n>min.recapt]
        n <- n[n>min.recapt]

        temp <- mean (maxd, na.rm = TRUE)

        if (!full) temp
        else {
            SE <- function(x) sqrt(var(x, na.rm=T)/sum(!is.na(x)))
            summaryd <- data.frame (Ncapt = names(table(n)),
                           n = as.numeric(table(n)),
                           mean = tapply(maxd, n, mean, na.rm=T),
                           se = tapply(maxd, n, SE))
            summaryd$mean[is.na(summaryd$mean)] <- NA  ## tidier
            summaryd <- rbind(summaryd, data.frame(Ncapt='Total', n=sum(!is.na(maxd)),
                mean=temp, se=SE(maxd)))
            list (MMDM = temp, data = data.frame(maxd = maxd, n=n), summary = summaryd )
        }
    }
}
############################################################################################

RPSV <- function (capthist, CC = FALSE)
{
    if (inherits (capthist, 'list')) {
        lapply(capthist, RPSV, CC)   ## recursive
    }
    else {
        RPSVx <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            x <- traps$x[cx]
            y <- traps$y[cx]
            n <- length(x)
            c(n = n-1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) - (sum(y))^2/n)
        }
        RPSVxy <- function (xy) {
            x <- xy[,1]
            y <- xy[,2]
            n <- length(x)
            c(n = n-1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) - (sum(y))^2/n)
        }
        if (nrow(capthist) < 1) return(NA)
        traps <- traps(capthist)
        det <- secr_expanddet(capthist)
        if (!all(det %in% .localstuff$individualdetectors))
            stop ("require individual detector type for RPSV")
        IDfactor <- getID(det, capthist)
        
        if (all(det %in% 'telemetry')) {
            lxy <- telemetryxy(capthist)
            if (is.null(lxy))
                temp <- NA
            else {
                temp <- lapply (lxy, RPSVxy)
            }
        }
        else if (all(det %in% .localstuff$polydetectors)) {
            
            if (is.null(xy(capthist)))
                temp <- NA
            else {
                lxy <- split ( xy(capthist), IDfactor)
                temp <- lapply (lxy, RPSVxy)
            }
        }
        else {
            if (any(det %in% 'telemetry')) {
                capthist <- subset(capthist, 
                                   occasions = det != 'telemetry', 
                                   traps = 1:(nrow(traps(capthist)-1)))
                IDfactor <- factor(animalID(capthist, names = TRUE), 
                    levels = rownames(capthist))
            }
            w <- split(trap(capthist, names = FALSE), IDfactor)
            temp <- lapply(w,RPSVx)
        }
        temp <- matrix(unlist(temp), nrow = 3)
        temp <- apply(temp,1,sum, na.rm=T)
        if (any(is.na(temp) | temp<0)) {
            temp <- NA   ## protected 2021-03-31
        }
        else {
            if (CC)
                temp <- sqrt((temp[2]+temp[3]) / (2 * temp[1]))
            else
                temp <- sqrt((temp[2]+temp[3]) / (temp[1]-1))
        }
        attr(temp,'names') <- NULL
        temp
    }
}

############################################################################################
## source ('d:\\density secr 1.3\\secr\\r\\mmdm.R')
## data(Peromyscus)

## MMDM(Peromyscus.WSG, full=T)$summary
##   Ncapt  n     mean        se
## 1     1  9       NA        NA
## 2     2  9 28.32839  9.434631
## 3     3 10 24.05921  9.335062
## 4     4  8 33.87949  5.471227
## 5     5  8 52.37655 15.470420
## 6     6  7 34.24929  5.961350
## 7 Total 42 33.93669  4.495646

## MMDM(Peromyscus.WSG, full=T)$summary[,3:4]/15.2
##       mean        se
## 1       NA        NA
## 2 1.863710 0.6206994
## 3 1.582843 0.6141488
## 4 2.228914 0.3599492
## 5 3.445826 1.0177908
## 6 2.253243 0.3921941
## 7 2.232677 0.2957662   <<<< 0.575?

## cf Otis et al 1978 p 87 Fig 23a

##################################################

ORL <- function (capthist, userdist = NULL, mask = NULL) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, ORL, userdist, mask)   
    }
    else {
        ORLx <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            as.numeric(max(distmat[cx, cx]))
        }
        ORLxy <- function (xy) {
            if (nrow(xy) == 1) 
                0
            else
                max(dist(cbind(xy$x, xy$y)))
        }
        if (nrow(capthist) < 1) return(NA)
        traps <- traps(capthist)
        det <- secr_expanddet(capthist)
        if (!all(det %in% .localstuff$individualdetectors))
            stop ("require individual detector type for ARL")
        distmat <- secr_valid.userdist(userdist, det, traps, traps, mask )
        prox  <- length(dim(capthist)) > 2
        IDfactor <- getID(det, capthist)
        
        if (all(det %in% 'telemetry')) {
            lxy <- telemetryxy(capthist)
            if (is.null(lxy))
                stop("no telemetry coordinates")
            else {
                maxd <- unlist(lapply (lxy, ORLxy))
                n <- unlist(lapply (lxy, nrow))
            }
        }
        else if (all(det %in% .localstuff$polydetectors)) {
            if (is.null(xy(capthist)))
                stop("no xy coordinates")
            else {
                lxy <- split (xy(capthist), IDfactor)
                maxd <- unlist(lapply (lxy, ORLxy))
                n <- unlist(lapply (lxy, nrow))
            }
        }
        else {
            if (any(det %in% 'telemetry')) {
                capthist <- subset(capthist, 
                    occasions = det != 'telemetry', 
                    traps = 1:(nrow(traps(capthist)-1)))
                IDfactor <- factor(animalID(capthist, names = TRUE), 
                    levels = rownames(capthist))
            }
            w <- split(trap(capthist, names = FALSE), IDfactor)
            maxd <- unlist(lapply(w, ORLx))
            n <- unlist(lapply(w, length))
        }
        data.frame(ORL = maxd, n = n)
    }
}
############################################################################################

## 2020-08-05, 2020-08-31

centroids <- function (capthist) {
    
    if (ms(capthist)) {
        nsess <- length(capthist)
        out <- lapply(capthist, centroids)
        cov <- lapply(capthist, covariates)
        ID <- unique(unlist(sapply(out, rownames)))
        IDxy <- array(
            dim = c(length(ID), 2, nsess), 
            dimnames = list(ID, c('meanx','meany'), session(capthist))
        )
        IDn <- matrix(0, 
            nrow = length(ID), 
            ncol = nsess,
            dimnames = list(ID, 1:nsess)
        )
        if (!is.null(cov[[1]])) {
            IDcov <- data.frame(
                matrix('', 
                    nrow = length(ID), 
                    ncol = ncol(cov[[1]]), 
                    dimnames = list(ID,names(cov[[1]]))
                )
            )
        }
            
        for (sess in 1:nsess) {
            IDxy[rownames(out[[sess]]), 1:2, sess] <- out[[sess]]
            IDn[rownames(out[[sess]]), sess] <- attr(out[[sess]], 'Ndetections')
        }
        attr(IDxy, 'Ndetections') <- IDn
        IDxy
    }
    else {
        centresx <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            x <- traps$x[cx]
            y <- traps$y[cx]
            c(n = length(x), meanx = mean(x, na.rm=TRUE), meany = mean(y, na.rm=TRUE))
        }
        centresxy <- function (xy) {
            x <- xy[,1]
            y <- xy[,2]
            c(n = length(x), meanx = mean(x, na.rm = TRUE), meany = mean(y, na.rm = TRUE))
        }
        if (nrow(capthist) < 1) return(NA)
        traps <- traps(capthist)
        det <- secr_expanddet(capthist)
        IDfactor <- getID(det, capthist)
        
        if (!all(det %in% .localstuff$individualdetectors))
            stop ("require individual detector type for centres")
        
        if (all(det %in% 'telemetry')) {
            lxy <- telemetryxy(capthist)  ## already list by animal
            if (is.null(lxy))
                temp <- NA
            else {
                temp <- lapply (lxy, centresxy)
            }
        }
        else if (all(det %in% .localstuff$polydetectors)) {
            
            if (is.null(xy(capthist)))
                temp <- NA
            else {
                lxy <- split ( xy(capthist), IDfactor)
                temp <- lapply (lxy, centresxy)
            }
        }
        else {
            if (any(det %in% 'telemetry')) {
                capthist <- subset(capthist, 
                    occasions = det != 'telemetry', 
                    traps = 1:(nrow(traps(capthist)-1)))
            }
            IDfactor <- factor(animalID(capthist, names = FALSE))
            w <- split(trap(capthist, names = FALSE), IDfactor)
            temp <- lapply(w,centresx)
        }
        out <- do.call(rbind, temp)
        rownames(out) <- rownames(capthist)
        ncapt <- out[,'n', drop = FALSE]
        out <- out[,c('meanx','meany')]
        attr(out, 'Ndetections') <- ncapt
        out
    }
}
############################################################################################

## 2025-06-16

t2r2 <- function (capthist)
{
    individualdetectors <- c('single','multi','proximity','count',
                             'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect',
                             'telemetry', 'capped')
    polydetectors <- c('polygon','transect','polygonX','transectX')
    getID <- function (det, capthist) {
        if (all(det %in% polydetectors)) {
            ID <- animalID(capthist, names = TRUE, sortorder = "ksn")
        }
        else {
            ID <- animalID(capthist, names = TRUE, sortorder = "snk")
        }
        factor(ID, levels = rownames(capthist))
    }
    
    if (inherits (capthist, 'list')) {
        lapply(capthist, t2r2)   ## recursive
    }
    else {
        squares <- function(x,y,n) {
            ssx <- sum(x^2) - (sum(x))^2/n
            ssy <- sum(y^2) - (sum(y))^2/n
            ssdx <- sum(diff(x)^2)
            ssdy <- sum(diff(y)^2)
            c(n = n-1, ssx = ssx, ssy = ssy, ssdx = ssdx, ssdy = ssdy)
        }        
        t2r2x <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            x <- traps$x[cx]
            y <- traps$y[cx]
            n <- length(x)
            squares(x,y,n)
        }
        t2r2xy <- function (xy) {
            x <- xy[,1]
            y <- xy[,2]
            n <- length(x)
            squares(x,y,n)
        }
        if (nrow(capthist) < 1) return(NA)
        traps <- traps(capthist)
        det <- secr_expanddet(capthist)
        if (!all(det %in% individualdetectors))
            stop ("require individual detector type for t2r2")
        IDfactor <- getID(det, capthist)
        
        if (all(det %in% 'telemetry')) {
            lxy <- telemetryxy(capthist)
            if (is.null(lxy))
                temp <- NA
            else {
                temp <- lapply (lxy, t2r2xy)
            }
        }
        else if (all(det %in% polydetectors)) {
            if (is.null(xy(capthist)))
                temp <- NA
            else {
                lxy <- split ( xy(capthist), IDfactor)
                temp <- lapply (lxy, t2r2xy)
            }
        }
        else {
            if (any(det %in% 'telemetry')) {
                capthist <- subset(capthist, 
                                   occasions = det != 'telemetry', 
                                   traps = 1:(nrow(traps(capthist)-1)))
                IDfactor <- factor(animalID(capthist, names = TRUE), 
                                   levels = rownames(capthist))
            }
            w <- split(trap(capthist, names = FALSE), IDfactor)
            temp <- lapply(w,t2r2x)
        }
        temp <- matrix(unlist(temp), nrow = 5)
        temp <- apply(temp,1,sum, na.rm=T)
        if (any(is.na(temp) | temp<0)) {
            temp <- NA   ## protected 2021-03-31
        }
        else {
            temp <- sum(temp[4:5])/sum(temp[2:3])
        }
        attr(temp,'names') <- NULL
        temp
    }
}
############################################################################################

## 2017-02-06 not exported
RPSVxy <- function (xy, CC = F) {
    x <- xy[,1]
    y <- xy[,2]
    n <- length(x)
    temp <- c(n = n-1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) - (sum(y))^2/n)
    if (CC)
        temp <- sqrt((temp[2]+temp[3]) / (2 * temp[1]))
    else
        temp <- sqrt((temp[2]+temp[3]) / (temp[1]-1))
    attr(temp,'names') <- NULL
    temp
}
##################################################

## not exported secr 4.4.5 2021-07-07

plotmoves <- function (capthist, byanimal = FALSE, withinsession = FALSE, 
    label = TRUE, arrows = TRUE, ...) {
    if (!ms(capthist)) {
        stop("plotmoves expects multi-session capthist")
    }
    cen <- centroids(capthist)
    ct <- attr(cen, 'Ndetections')>0
    ok <- apply(ct,1,sum)>1   # at least 2 sessions
    if (!byanimal) {
        plot(traps(capthist[[1]]), ...)  
    }
    for (j in which(ok)) {
        ch <- suppressWarnings(subset(capthist, rownames(cen)[j]))
        if (byanimal) {
            plot(traps(capthist[[1]]), ...)
            mtext(side=3, rownames(cen)[j], line=0.4, cex=0.7)
        }
        if (withinsession) {
            plot(ch, add=T, tracks = T, varycol=FALSE, title='', subtitle='')
        }
        for (i in 1:4) {
            move2 <- (cen[j,1,i]-cen[j,1,i+1])^2 + (cen[j,2,i+1]-cen[j,2,i])^2 
            if (!is.na(move2) && move2>0.001 && arrows) {
                arrows(cen[j,1,i], cen[j,2,i], cen[j,1,i+1], cen[j,2,i+1], 
                    lwd = 1.5, angle = 15, length=0.15)
            }
            else {
                segments(cen[j,1,i], cen[j,2,i], cen[j,1,i+1], cen[j,2,i+1], 
                    lwd = 1.5)
            }
        }
        if (label) {
            points(cen[j,1,], cen[j,2,], pch = 16, col='yellow', cex=2)
            text(cen[j,1,], cen[j,2,], 1:5, cex=0.9)
        }
    }
    d <- apply(cen, 1, function (xy) (diff(xy[1,])^2 + diff(xy[2,])^2)^0.5)
    d[!is.na(d)]  # vector of consecutive moves
}

# par(mfrow=c(3,8), mar=c(2,2,2,2))
# d <- plotmoves(ovenCHp, label=T, arrows=F)
# symbols(circles=80, )
