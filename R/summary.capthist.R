#############################################################################
## package 'secr'
## summary.capthist.R
## 2019-04-04 moved from methods.R 
## 2019-11-10 summary of covariates
#############################################################################

summary.capthist <- function(object, terse = FALSE, moves = FALSE, tpa = FALSE, ...) {
    ## recursive if list of capthist
    if (ms(object)) {
        if (terse) {
            sapply (object, summary, terse = TRUE, moves = moves, tpa = tpa, ...)
        }
        else {
            lapply (object, summary, terse = FALSE, moves = moves, tpa = tpa, ...)
        }
    }
    else {
        object <- check3D(object)
        trps <- traps(object)
        nd <- ndetector(trps)
        if (terse) {   ## 2017-11-06, 2019-01-22
            c(Occasions = ncol(object),
              Detections = sum(abs(object)),
              Animals = nrow(object),
              Detectors = nd,
              Moves = if (moves) sum(unlist(sapply(moves(object), function(y) y>0))) else NULL,
              Animals2 = if (tpa) sum(trapsPerAnimal(object)[-1]) else NULL
            )
        }
        else {
            detector <- expanddet(object) # detector(traps)
            cutval <- attr(object, 'cutval',exact = TRUE)   # signal strength only
            
            # ni, ui, fi, M, losses etc.
            nocc <- ncol(object)
            counts <- matrix(0, nrow = 8, ncol = nocc)
            signalsummary <- NULL
            if (nrow(object) > 0) {
                tempx <- apply( object[,,,drop=F], c(1,2), function(x) sum(abs(x))>0)
                tempx3 <-  apply( object[,,,drop=F], c(1,2), function(x) any(x<0))
                if (nocc>1) {  # distinction may not be needed...
                    counts [1,] <- apply(tempx, 2, function(x) sum(abs(x)>0) )
                    tempx2 <- apply(tempx, 1, function(x) cumsum(abs(x))>0)
                    counts [4,] <- apply(tempx2,1,sum)
                    counts [2,] <- c(counts[4,1],diff(counts[4,]))
                    counts [3,] <- tabulate(apply(tempx,1, function(x) sum(abs(x)>0)),nbins = nocc)
                    counts [5,] <- apply(tempx3, 2, sum)
                }
                else {
                    counts [1,1] <- sum(abs(tempx)>0)
                    counts [2,1] <- counts [1,1]
                    counts [3,1] <- counts [1,1]
                    counts [4,1] <- counts [1,1]
                    counts [5,1] <- sum(tempx<0)
                }
                counts [6,] <- apply(object[,,, drop=F], 2, function(x) sum(abs(x)))
                tempt <- apply(object[,,,drop=F],c(2,3), function(x) sum(abs(x))>0)
                counts [7,] <- apply(tempt,1,sum)
            }
            if (!is.null(trps)) {
                if (!is.null(usage(trps))) {
                    counts[8,] <- apply(usage(trps),2,function(x) sum(x>0))
                }
                else {
                    counts[8,] <- rep(nd, nocc)
                }
            }
            counts <- as.data.frame(counts)
            dimnames(counts) <- list(c('n','u','f','M(t+1)','losses','detections',
                                       'detectors visited','detectors used'), 1:nocc)
            counts$Total <- apply(counts, 1, sum)
            counts$Total[4] <- counts[4, nocc]
            
            PSV <-  NULL
            dbar <- NULL
            movesummary <- NULL
            tpasummary = NULL
            
            if (is.null(trps)) {
                trapsum <- NULL
                signalsummary <- NULL
            }
            else {
                
                if (all(detector %in% .localstuff$individualdetectors)) {
                    if (-diff(counts$Total[1:2]) > 1)
                        PSV <- RPSV(object)
                    if (all(detector(trps) %in% .localstuff$exclusivedetectors))
                        dbar <- dbar(object)
                    if (moves) {
                        mov <- moves(object)
                        mov <- lapply(mov, function(x) x[x>0])
                        movesummary <- list (peranimal = table(sapply(mov, length)),
                                             distance = summary(unlist(mov)))
                    }
                    if (tpa) {
                        tpasummary <- trapsPerAnimal(object)
                    }
                }
                trapsum <- summary(trps)
                if (all(detector == 'signal'))
                    signalsummary <- summary(signal(object))
                if (all(detector == 'signalnoise'))
                    signalsummary <- list(signal = summary(signal(object)),
                                          noise = summary(noise(object)),
                                          diffSN = summary(signal(object)-noise(object)))
            }
            zeros <- sum(apply(abs(object),1,sum)==0)
            
            xyl <- telemetryxy(object)
            if (is.null(xyl))
                telemsummary <- NULL
            else {
                telemocc <- which(detector(trps) == 'telemetry')
                zeros <-  sum(apply(abs(object)[,-telemocc,, drop = FALSE],1,sum)==0)
                ntelem <- sapply(xyl, nrow)
                if (nrow(object) == 0 | ncol(object)==1 | all(detector(trps) == 'telemetry'))
                    nteldet <- 0
                else {
                    s1 <- apply(abs(object)[,-telemocc,, drop = FALSE]>0,1,any)
                    nteldet <- sum(s1 [row.names(object) %in% names(xyl)])
                }
                telemsummary <- c(n=length(xyl), ndet=nteldet, min=min(ntelem),
                                  max=max(ntelem), mean=mean(ntelem), sd = sd(ntelem))
                counts[7:8, telemocc] <- 0  ## no detector visits/used
                counts$Total[7:8] <- counts$Total[7:8]-length(telemocc)
            }
            
            ## resightings
            markocc <- markocc(trps)
            sighting <- sighting(trps)
            if (sighting) {
                sightings <- matrix(0, nrow = 5, ncol = nocc+1,
                                    dimnames = list(c('ID','Not ID','Unmarked','Uncertain',
                                                      'Total'), c(1:nocc, 'Total')))
                Tm <- Tm(object)
                Tu <- Tu(object)
                Tn <- Tn(object)
                unresolved <- markocc == -1
                sightings[1, c(markocc < 1, FALSE)] <- unlist(counts[6,c(markocc < 1, FALSE)])
                if (!is.null(Tm)) {
                    if (is.matrix(Tm))
                        sightings[2, 1:nocc] <- apply(Tm,2,sum)
                    else {
                        sightings[2, (1:nocc)[markocc<1]] <- NA
                        sightings[2, nocc+1] <- Tm
                    }
                }
                if (!is.null(Tu)) {
                    if (is.matrix(Tu))
                        sightings[3, 1:nocc] <- apply(Tu,2,sum)
                    else {
                        sightings[3, (1:nocc)[markocc<1]] <- NA
                        sightings[3, nocc+1] <- Tu
                    }
                }
                if (!is.null(Tn)) {
                    if (is.matrix(Tn))
                        sightings[4, 1:nocc] <- apply(Tn,2,sum, na.rm=T)
                    else {
                        sightings[4, (1:nocc)[markocc<1]] <- NA
                        sightings[4, nocc+1] <- Tn
                    }
                }
                ## avoid replacing summed counts
                sightings[, nocc + 1] <- ifelse (sightings[, nocc + 1]==0,
                                                 apply(sightings, 1, sum),
                                                 sightings[, nocc + 1])
                
                sightings[5, 1:(nocc+1)] <- apply(sightings[1:4, 1:(nocc+1), drop = FALSE], 2, sum, na.rm=T)
            }
            else sightings <- NULL
            
            if (!is.null(attr(object, 'nontarget',exact = TRUE))) {
                nontarg <- attr(object, 'nontarget',exact = TRUE)
                nontarget <- data.frame(as.list(apply(nontarg, 2, sum)), sum(nontarg))
                names(nontarget) <- c(1:nocc, 'Total')
                rownames(nontarget) <- 'detectors nontarget'
            }
            else nontarget <- NULL
            
            if (is.null(trps)) {
                counts <- counts[1:6,]
            }
            
            tempcovar <- covariates(object)
            if (!is.null(tempcovar) && ((nrow(tempcovar)>0) & (ncol(tempcovar)>0))) {
                covsummary <- summary(stringsAsFactors(tempcovar))   ## force character to factor 2020-07-14 
            }
            else {
                covsummary <- NULL
            }
            temp <- list (
                detector = detector,
                ndetector = nd,
                trapsum = trapsum,
                counts = counts,
                zeros = zeros,
                dbar = dbar,
                RPSV = PSV,
                moves = movesummary,
                tpa = tpasummary,
                cutval = cutval,        # signal, signalnoise only
                signalsummary = signalsummary,
                telemsummary = telemsummary,
                sightings = sightings,
                nontarget = nontarget,
                covsummary = covsummary
            )
            class(temp) <- 'summary.capthist'
            temp
        }
    }
}
############################################################################################

counts <- function (CHlist, counts = 'M(t+1)') {
    if (!inherits(CHlist, 'capthist'))
        stop ("require capthist object")
    getc <- function (cnt) {
        getcnt <- function(x, maxocc) {
            temp <- x$counts[cnt,]
            lt <- length(temp)
            matrix(c(temp[-lt], rep(NA, maxocc-lt+1), temp[lt]), nrow = 1)
        }
        if (!is.list(CHlist))
            summary(CHlist)$counts[cnt,]
        else {
            maxocc <- max(sapply(CHlist,ncol))
            abind(lapply(summary(CHlist), getcnt, maxocc), along=1,
                  new.names=list(session(CHlist), c(1:maxocc, 'Total')))
        }
    }
    temp <- lapply (counts, getc)
    names(temp) <- counts
    temp
}

print.summary.capthist <- function (x, ...) {
    cat ('Object class      ', 'capthist', '\n')
    nonspatial <- is.null(x$trapsum)
    if (!nonspatial)
        print(x$trapsum, terse=TRUE)
    
    cat ('\nCounts by occasion \n')
    print(x$counts, ...)
    if (x$zeros>0)
        cat ('\nEmpty histories : ', x$zeros, '\n')
    
    if (!is.null(x$moves)) {
        cat ('\nNumber of movements per animal')
        print(x$moves$peranimal)
        cat ('\nDistance moved, excluding zero (m)\n')
        print(x$moves$distance)
    }
    
    if (!is.null(x$tpa)) {
        cat ('\nNumber of detectors per animal\n')
        print(x$tpa)
    }
    
    if (!nonspatial) {
        if (all(x$detector %in% c('signal', 'signalnoise'))) {
            cat ('Signal threshold ', x$cutval, '\n')
            print (x$signalsummary)
        }
        if (!is.null(x$telemsummary)) {
            cat (x$telemsummary[1], "telemetered animals,", x$telemsummary[2], "detected\n")
            cat (paste(x$telemsummary[3:4], collapse="-"), "locations per animal, mean = ",
                 paste(round(x$telemsummary[5:6],2), collapse=", sd = "), "\n")
        }
    }
    if (!is.null(x$sightings)) {
        cat ('\nSightings by occasion \n')
        print(x$sightings, ...)
        cat ('\n')
    }
    if (!is.null(x$nontarget)) {
        cat ('\nNon-target interference by occasion \n')
        print(x$nontarget, ...)
        cat ('\n')
    }
    if (!is.null(x$covsummary)) {
        cat ('\nIndividual covariates\n')
        print (x$covsummary)
        cat('\n')
    }
    
}
############################################################################################
