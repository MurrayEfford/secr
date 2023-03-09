## package 'secr'
## sim.capthist.R
## simulate capture histories
###############################################################################

expands <- function (x, s, default = 1) {
    if (is.null(x))
        rep(default, s)
    else {
        y <- numeric(s)
        y[] <- x
        y
    }
}

expandsk <- function (x, s, k, default = 1) {
    if (is.null(x))
        matrix(default, nrow = s, ncol = k)
    else
        matrix(x, nrow = s, ncol = k)
}

sim.capthist <- function (
    traps,
    popn = list(D = 5, buffer = 100, Ndist = 'poisson'),
    detectfn = 0,
    detectpar = list(),
    noccasions = 5,
    nsessions = 1,
    binomN = NULL,
    exactN = NULL,
    p.available = 1,
    renumber = TRUE,
    seed = NULL,
    maxperpoly = 100,
    chulltol = 0.001,
    userdist = NULL,
    savepopn = FALSE
    )

#
# Simulate detection of known population
#

## A note on sort order of records  2009 12 01

## Simulation routines return the primary detection data in 'ski' order,
## meaning that occasion (s) changes fastest and individual (i) slowest -
## this allows efficient sequential output as new animals are detected.
## In R the order is changed to 'isk' as this is pictorially more natural.

## Secondary data (xy locations, signal strength) are generated in the order
## 'kis' (detector (k) changing fastest), because the simulation routines use -
## for (s=0; s<*ss; s++)
##   for (i=0; i<*N; i++)
##     for (k=0; k<*kk; k++) {
##     ...
##     }
## or similar loops.
## For consistency, all data are sorted to 'isk' order before returning to R.
## Secondary simulated data are re-sorted into 'ksi' order by creating the
## index 'start' within C as required for prwipolygon & prwitransect.
## The functions 'animalID', 'occasion', and 'trap' are the safest way to
## retrieve values for detections in isk order.

{
    poplist <- inherits(popn,'popn') & inherits(popn,'list')

    # moved here 2020-12-06
    if (is.character(detectfn)) {
        detectfn <- detectionfunctionnumber(detectfn)
    }
    
    ##------------------------------------------------------------------------
    ## multi-session loop
    if (poplist | (nsessions > 1)) {
        
        if (poplist & (nsessions>1) & (length(popn) != nsessions))
            stop ("incompatible use of popn list and nsessions>1")

        nsessions <- ifelse (poplist, length(popn), nsessions)
        if (ms(traps) & (length(traps) != nsessions))
            stop ("incompatible use of traps list and nsessions>1")

        output <- vector(nsessions, mode='list')
        nocc <- numeric(nsessions)
        nocc[] <- noccasions

        ##################
        ## set random seed
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
        #####################################################################
        ## generate population list if not provided
        if (!inherits(popn,'popn')) {
            popn <- replacedefaults(list(D = 5, buffer = 100,
                Ndist = 'poisson'), popn)
            ## will fail with multiple traps
            popn <- sim.popn (popn$D, core = traps, buffer = popn$buffer,
                covariates = NULL, Ndist = popn$Ndist, nsessions = nsessions)
            poplist <- TRUE
        }

        #####################################################################
        ## availability preliminaries
        if (poplist) {
            if (any(p.available != 1))
                warning ("incomplete availability not implemented ",
                         "for population lists")
            if(renumber)
                warning ("typically use renumber = FALSE for multiple sessions")
        }
        else {
            if (!(length(p.available) %in% 1:2))
                stop ("p.available must be vector of 1 or 2 probabilities")
            availability <- 'random'
            if (length(p.available) == 1)
                ## random temporary emigration
                available <- runif(nrow(popn)) < p.available
            else {
                ## Markovian temporary emigration
                availability <- 'Markov'
                equilibrium.p <- (p.available[2] / (1-p.available[1]+p.available[2]))
                available <- runif(nrow(popn)) < equilibrium.p
            }
        }
        #####################################################################
        ## session-specific detection 2015-04-01; extended to HPX 2021-03-25
        if (detectfn %in% c(0:7, 14:20)) {
            df0name <- if (detectfn %in% (0:7)) 'g0' else 'lambda0'
            dfzname <- if (detectfn %in% (5:6)) 'x' else 'z'
            df0 <- expands(detectpar[[df0name]], nsessions, default = NULL)
            sigma <- expands(detectpar[['sigma']], nsessions, default = NULL)
            dfz <- expands(detectpar[[dfzname]], nsessions, default = NULL)
        }
        # detectfn, 9, c('g0', 'sigma'), c('b0', 'b1'))
        # detectfn, 10:13, c('g0', 'sigma'), c('beta0','beta1'))
        
        ########################################################################
        ## loop over sessions
        for (t in 1:nsessions) {
            ## if (R > 1)  ## modified 2015-10-29

            if (poplist)
                temppop <- popn[[t]]
            else {
                temppop <- subset(popn, available)
                ##-------------------------------------------------------------
                ## update availability in preparation for next session
                if (availability == 'random') {
                    available <- runif(nrow(popn)) < p.available
                }
                else {
                    p.vect <- ifelse(available, p.available[1], p.available[2])
                    available <- runif(nrow(popn)) < p.vect
                }
                #--------------------------------------------------------------
            }

            ## select session-specific traps if necessary
            trps <- if (ms(traps)) traps[[t]] else traps

            ## session-specific detection parameters; extended to HPX 2021-03-25
            if (detectfn %in% c(0:7, 14:20)) {
                detectpar[[df0name]] <- df0[t]
                detectpar[['sigma']] <- sigma[t]
                detectpar[[dfzname]] <- dfz[t]
            }
            # detectfn, 9, c('g0', 'sigma'), c('b0', 'b1'))
            # detectfn, 10:13, c('g0', 'sigma'), c('beta0','beta1'))

            ## recursive call for each session
            ## we have set the seed above
            output[[t]] <- sim.capthist(trps, temppop, detectfn, detectpar,
                noccasions = nocc[t], nsessions = 1,
                binomN = binomN, exactN = exactN, p.available = 1,
                renumber = renumber, seed = NULL, maxperpoly)

            session( output[[t]] ) <- t   ## added 2011-09-09
        }
        ########################################################################

        class(output) <- c('capthist', 'list')
        names(output) <- 1:nsessions
        output
    }   ## end of multi-session

    else {
        ##################
        ## set random seed
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
        ##################
        if (is.null(detector(traps)))
            stop ("'traps' lacks detector type")
        if (length(detector(traps)) != noccasions) 
            detector(traps) <- rep(detector(traps)[1], noccasions)
        
        usge <- usage(traps)
        if (is.null(usge))
            usge <- matrix (1, nrow = ndetector(traps), ncol = noccasions)
        else {
            if (nrow(usge) != ndetector(traps))
                stop ("invalid usage matrix; number of rows ",
                      "must match number of detectors")
            if (ncol(usge) != noccasions) {
                noccasions <- ncol(usge)
                warning ("'noccasions' does not match usage ",
                         "attribute of 'traps'; ignored")
            }
        }

        validatepar <- function (x, xrange) {
            xname <- deparse(substitute(x))
            if (is.null(x))
                stop ("no value for ", xname)
            if (any(is.na(x)))
                stop ("NA is not a valid value for ", xname)
            if (any(x < xrange[1]))
                warning ("value for ", xname, " is less than minimum ",
                    xrange[1])
            if (any(x > xrange[2]))
                warning ("value for ", xname, " is greater than maximum ",
                    xrange[2])
        }

        if (all(detector(traps) %in% c('signal'))) {
            if (!(detectfn %in% c(10,11))) {
                warning ("forcing detection function = 10 for signal detectors")
                detectfn <- 10
            }
        }
        else if (all(detector(traps) %in% c('signalnoise'))) {
            if (!(detectfn %in% c(12,13))) {
                warning ("forcing detection function = 12 for signalnoise detectors")
                detectfn <- 12
            }
        }
        
        if (any(detector(traps) %in% .localstuff$polydetectors)) {
           if (!all(detectfn %in% 14:20)) {
            stop("polygon and transect detectors use hazard detection functions 14:20 (HHN etc.)")
           }
        }

        ## Detection function parameters

        ##    0  halfnormal
        ##    1  hazard rate
        ##    2  exponential
        ##    3  compound halfnormal
        ##    4  uniform
        ##    5  w-exponential
        ##    6  annular normal
        ##    7  cumulative lognormal
        ##    8  cumulative gamma??
        ##    9  binary signal strength (b0 = (beta0-c)/sdS, b1 = beta1/sdS)
        ##   10  signal strength (signal detectors only)
        ##   11  signal strength with spherical spreading (signal detectors only)
        ##   12  signal-noise (signalnoise detectors only)
        ##   13  signal-noise with spherical spreading (signalnoise detectors only)
        ##   14  hazard halfnormal
        ##   15  hazard hazard-rate
        ##   16  hazard exponential
        ##   17  hazard annular normal
        ##   18  hazard cumulative gamma
        ##   19  hazard variable power
        ##   20  hazard pixelar 2021-03-25

        if (!is.null(usage(traps))) {
            if (!is.null(noccasions)) {
                if (noccasions != ncol(usage(traps)))
                    warning ("noccasions differs from ncol(usage); using latter")
            }
            noccasions <- ncol(usage(traps))
        }
        if (detectfn %in% c(0:4))  defaults <- list(g0 = 0.2, sigma = 25, z = 1)
        if (detectfn %in% c(5))    defaults <- list(g0 = 0.2, sigma = 25, w = 10)
        if (detectfn %in% c(6))    defaults <- list(g0 = 0.2, sigma = 25, w = 10)
        if (detectfn %in% c(7,8))    defaults <- list(g0 = 0.2, sigma = 25, z = 5)
        if (detectfn %in% c(9))  defaults <- list(b0 = 1, b1=-0.1, cutval = 60,
            tx = 'identity')
        if (detectfn %in% c(10,11))  defaults <- list(beta0 = 90, beta1=-0.2,
            sdS = 2, cutval = 60, sdM = 0, tx = 'identity')
        if (detectfn %in% c(14:16,19))  defaults <- list(lambda0 = 0.2, sigma = 25, z = 1)
        if (detectfn %in% c(17))  defaults <- list(lambda0 = 0.2, sigma = 25, w =10)
        if (detectfn %in% c(18))  defaults <- list(lambda0 = 0.2, sigma = 25, z = 5)
        
        # 2021-03-25
        if (detectfn %in% c(20))  defaults <- list(lambda0 = 0.2, sigma = spacing(traps)/2, z = 1)
        
        if (detectfn %in% c(12,13))  defaults <- list(beta0 = 90, beta1=-0.2,
            sdS = 2, cutval = 10, muN = 40, sdN = 2, sdM = 0, tx = 'identity')
        else defaults <- c(defaults, list(truncate = 1e+10, recapfactor = 1.0))

        if (is.null(binomN)) {
            detectpar$binomN <- 0  ## default Poisson counts
        }
        else {
            detectpar$binomN <- if (tolower(binomN[1]) == 'usage') 
                1       ## code for 'binomial size from usage' 2012-12-22
            else 
                binomN  ## as input
        }
        
        detectpar <- replacedefaults(defaults, detectpar)

        # extended to HPX 2021-03-25
        if (detectfn %in% c(0:7, 14:20)) {
            if (detectfn %in% (0:7)) {
                g0    <- expandsk(detectpar$g0, s = noccasions, k = ndetector(traps))
                if (any(detector(traps) %in% .localstuff$countdetectors)) {
                    if (detectpar$binomN == 0)
                        validatepar(g0, c(0,Inf))  ## Poisson lambda
                    else
                        validatepar(g0, c(0,1))    ## Binomial p
                }
                df0 <- g0
            }
            else  {
                lambda0 <- expandsk(detectpar$lambda0, s = noccasions, k = ndetector(traps))
                validatepar(lambda0, c(0,Inf))  ## Poisson lambda
                df0 <- lambda0
            }

            sigma <- matrix(expands(detectpar$sigma, noccasions), ncol = 1)
            z     <- expands(ifelse(detectfn %in% c(5,6),
                detectpar$w, detectpar$z), noccasions)
            validatepar(sigma[,1], c(1e-10,Inf))
            validatepar(z, c(0,Inf))
        }
        # detectfn, 9, c('g0', 'sigma'), c('b0', 'b1'))
        # detectfn, 10:13, c('g0', 'sigma'), c('beta0','beta1'))

        # Acoustic detection function parameters
        if (detectfn %in% c(10,11,12,13)) {
            df0 <- 1   ## dummy value - not used
            tx <- detectpar$tx
            cutval <- detectpar$cutval
            sdM <- detectpar$sdM
            beta0 <- expands(detectpar$beta0, noccasions)
            beta1 <- expands(detectpar$beta1, noccasions)
            sdS   <- expands(detectpar$sdS, noccasions)
            muN   <- expands(detectpar$muN, noccasions)
            sdN   <- expands(detectpar$sdN, noccasions)
            validatepar(beta0, c(-Inf,Inf))
            validatepar(beta1, c(-Inf,Inf))
            validatepar(sdS, c(0,Inf))
            validatepar(sdN, c(0,Inf))
        }
        else if (detectfn %in% c(9)) {
            df0 <- expandsk(detectpar$b0, s = noccasions, k = ndetector(traps))
            sigma <- expands(detectpar$b1, noccasions)
            z <- 0
            cutval <- detectpar$cutval
        }
        else {
            cutval <- NULL
            truncate <- ifelse(is.null(detectpar$truncate),
                1e+10, detectpar$truncate)
            validatepar(truncate, c(1e-10, Inf)) ## must be positive
        }

        # Allow general learned response for traps
        if (any(detector(traps) %in% c('single','multi'))) 
        {
            rfl <- length(detectpar$recapfactor)
            if (!(rfl %in% 1:2))
                stop ("invalid recapfactor")
            if (rfl==1)
                detectpar$recapfactor <- c(detectpar$recapfactor,1)
            df0 <- cbind(df0, df0 * detectpar$recapfactor[1])
            sigma <- cbind(sigma, sigma * detectpar$recapfactor[2])
        }
     
        else {
            if (any(detectpar$recapfactor != 1.0))
                stop("learned response available only for 'single', 'multi' detectors")
        }

        #-----------------------------------------------------------------------------#
        if (!inherits(popn,'popn')) # generate if not provided
        {
            popn <- replacedefaults(list(D = 5, buffer = 100,
                Ndist = 'poisson'), popn)
            popn <- sim.popn (popn$D, core = traps, buffer = popn$buffer,
                covariates = NULL, Ndist = popn$Ndist)
        }
        #-----------------------------------------------------------------------------#

        ## user-provided distances
        if (is.null(userdist)) {
            distmat2 <- getdistmat2(traps, popn, NULL, detectfn==20)
        }
        else {
            ## move towards IHP/k simulations
            ## mask is NULL unless IHP or linear
            mask <- attr(popn, 'mask')
            ## ASSUME MASK HAS NONEUC REQUIREMENTS 2014-10-30
            ## fails if any(detector %in% .localstuff$polydetectors)
            distmat2 <- valid.userdist(userdist,
                                  detector(traps),
                                  xy1 = traps,
                                  xy2 = popn,         # animals, 2017-04-06
                                  mask = mask)^2
        }
        #-----------------------------------------------------------------
  
        stopiferror <- function(resultcode, fnname) {
            if (resultcode[1] != 0)
                stop ("call to '", fnname[1], "' failed")
        }
        ##-----------------------------------------------------------------
        
        N <- nrow(popn)
        animals <- as.matrix(popn)
        #k <- nrow(traps)
        k <- ndetector(traps)
        dettype <- detectorcode(traps, noccasions = noccasions)
        HPXpoly <- (detectfn == 20) && any(detector(traps) %in% .localstuff$polydetectors)
        if (HPXpoly)
            simfunctionname <- 'trappingproximity'
        else
            simfunctionname <- paste0('trapping', detector(traps))
        if (length(simfunctionname)==1 & noccasions>1)
            simfunctionname <- rep(simfunctionname, noccasions)
        if (length(simfunctionname)>1 & length(simfunctionname)!=noccasions)
            stop("provide one detector name for each occasion")

        ## vector binomN
        detectpar$binomN <- expandbinomN(detectpar$binomN, dettype)
    
        
        if (detectfn %in% c(0:7, 14:20)) {
            trappingargs <- list(
            simfunctionname[1],                 # 1
            as.double(df0),                     # 2
            as.double(sigma),                   # 3
            as.double(z),                       # 4
            as.matrix(distmat2),                # 5   assume ntraps x nanimals
            as.matrix(usge),                    # 6
            as.integer(detectfn),               # 7 
            as.double(truncate^2),              # 8
            as.integer(detectpar$binomN)        # 9 
            )
        }
        ## 'count' includes presence
        ##-----------------------------------------------------------------
        ## remember: output from trappingxxxx is sorted so captured animals
        ## precede all others
        ## caught[] contains capture sequence for each of N animals
        
        if ((all(detector(traps) %in% c('single','multi','proximity','count', 'capped')))
            || detectfn==20) {
            
            ## BUG IN THIS CODE 2016-10-18 Observations assigned to wrong detectors
            # if (length(unique(detector(traps)))==1) {
            #     sm <- detector(traps)[1] %in% c('single','multi')
            #     k1 <- if (sm) 1 else k
            #     trappingargs[[18]] <- integer (N * noccasions * k1)  ## output
            #     temp <- do.call(.C, trappingargs)
            #     stopiferror(temp$resultcode, simfunctionname)
            #     if (temp$n > 0) {
            #         if (sm) {
            #             temp$value[temp$value==0] <- NA
            #             occ <- rep(0:(noccasions-1), temp$n)
            #             trp <- temp$value[1:(temp$n * noccasions)]-1
            #             id <- rep(which(temp$caught>0)-1, each = noccasions)
            #             index <- (id * k + trp) * noccasions + occ + 1
            #             w[index] <- 1
            #         }
            #         else {
            #             occ <- rep(0:(noccasions-1), temp$n*k)
            #             trp <- rep(rep(0:(k-1), each = noccasions), temp$n)
            #             id <- rep(which(temp$caught>0)-1, each = noccasions*k)
            #             index <- (id * k + trp) * noccasions + occ + 1
            #             w[index] <- temp$value[1:(temp$n * k * noccasions)]
            #         }
            #     }
            # }
            # else 
            w <- array(0, dim=c(N, noccasions, k), dimnames = list(rownames(popn), 1:noccasions, 1:k))
            
            ## 2017-12-02
            ## external tracking of previous captures is needed to allow single-occasion calls to trappingxxx
            ## functions; I hope this also works with multi-occasion calls (currently disabled)
            lastcapt <- rep(0, N)
            for(s in 1:noccasions)   ## about 25% slower with one occasion at a time
            {
                sm <- dettype[s] %in% c(-1,0)    # single, multi
                k1 <- if (sm) 1 else k
                trappingargs[[1]] <- simfunctionname[s]
                trappingargs[[2]] <- as.double(df0[s,])
                trappingargs[[3]] <- as.double(sigma[s,])
                trappingargs[[6]] <- as.matrix(usge[,s, drop = FALSE])    ## usage on occasion s
                trappingargs[[9]] <- as.integer(detectpar$binomN[s])
                temp <- do.call(simfunctionname[s], trappingargs[-1])
                stopiferror(temp$resultcode, simfunctionname[s])
                # if any animals caught...
                if (temp$n > 0) {
                    if (sm) {
                        # maximum 1 capture/animal/occasion
                        occ <- rep(s, sum(temp$caught>0))
                        trp <- temp$value[temp$caught]
                        id <- (1:N)[temp$caught>0]
                        w[cbind(id, occ, trp)] <- 1
                        lastcapt[temp$caught>0] <- s    ## for recap
                    }
                    else {
                        occ <- rep(s, N*k)
                        trp <- rep(1:k, each = N)
                        id <- rep(1:N, k)
                        w[cbind(id, occ, trp)] <- temp$value
                        # lastcapt not updated because recapfactor does not apply
                    }
                                    }
            }
            ## drop empty histories
            w <- w[apply(w,1,sum)>0,,, drop = FALSE]
            class(w) <- 'capthist'
            traps(w) <- traps

            if (temp$n > 0 && HPXpoly) {
                ## put XY coordinates in attribute
              ## xy(w) <- data.frame(animals[animalID(w, names = TRUE),])
              xy(w) <- data.frame(animals[animalID(w, names = TRUE, sortorder = 'ksn'),])
            }
            else {
                xy(w) <- NULL
            }
            
            
        }
        ##-----------------------------------------------------------------------
        else if (detector(traps)[1] %in% c('polygonX','transectX')) {
            if (detector(traps)[1] == 'polygonX') {
                polynames <- levels(polyID(traps))
                nk <- length(polynames)
                k <- table(polyID(traps))
                temp <- trappingpolygonX (
                    as.double(df0),
                    as.double(sigma),
                    as.double(z),
                    as.integer(nk),
                    as.integer(k),
                    as.matrix(animals),
                    as.matrix(traps),
                    as.matrix(usge),
                    as.integer(detectfn),
                    as.double(truncate^2),
                    as.integer(detectpar$binomN)
                )
            }
            else {
                polynames <- levels(transectID(traps))
                nk <- length(polynames)
                if (nk>1) stop("simulation not working for multiple transects in secr 4.5")
                k <- table(transectID(traps))
                temp <- trappingtransectX (
                    as.double(df0),
                    as.double(sigma),
                    as.double(z),
                    as.integer(nk),
                    as.integer(k),
                    as.matrix(animals),
                    as.matrix(traps),
                    as.matrix(usge),
                    as.integer(detectfn),
                    as.double(truncate^2)
                )
            }
            if (temp$resultcode != 0)
                stop ("call to '", simfunctionname, "' failed")
            w <- array(0, dim=c(N, noccasions, nk), 
                       dimnames = list(rownames(popn), 1:noccasions, polynames))

            if (temp$n > 0) {
                # maximum 1 capture/animal/occasion
                # sorted by animal, occasion
                occ <- rep(1:noccasions, temp$n)
                trp <- temp$value[1:(temp$n*noccasions)]
                id <- rep(1:temp$n, each = noccasions)
                w[cbind(id, occ, trp)] <- 1
            }
            ## drop empty histories
            w <- w[apply(w,1,sum)>0,,, drop = FALSE]

            class(w) <- 'capthist'  
            traps(w) <- traps
            
            if (temp$n > 0) {
                ## put XY coordinates in attribute
                nd <- sum(abs(w) > 0)
                detectedXY <- data.frame(matrix(ncol = 2,
                    temp$detectedXY[1:(2*nd)]))
                names(detectedXY) <- c('x','y')
                xy(w) <- detectedXY
            }
            else
                xy(w) <- NULL
            }
        ##-----------------------------------------------------------------------
        else if (detector(traps)[1] %in% c('polygon','transect','telemetry')) {
            if (detector(traps)[1] == 'polygon') {
                nk <- length(levels(polyID(traps)))
                k <- table(polyID(traps))
                
                temp <- trappingpolygon (
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(nk),
                           as.integer(k),
                           as.matrix(animals),
                           as.matrix(traps),
                           as.matrix(usge),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           as.integer(detectpar$binomN),
                           as.integer(maxperpoly)
                           )
                polynames <- levels(polyID(traps))
            }
            ##-----------------------------------------------------------------------
            if (detector(traps)[1] == 'transect') {
                nk <- length(levels(transectID(traps)))
                if (nk>1) stop("simulation not working for multiple transects in secr 4.5")
                k <- table(transectID(traps))
                temp <- trappingtransect (
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.integer(nk),
                           as.integer(k),
                           as.matrix(animals),
                           as.matrix(traps),
                           as.matrix(usge),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           as.integer(detectpar$binomN),
                           as.integer(maxperpoly)
                           )
                polynames <- levels(polyID(traps))
            }
            ##-----------------------------------------------------------------------
            if (detector(traps)[1] == 'telemetry') { 
                nk <- 1
                polynames <- '1'
                maxperpoly <- 2000
                if (is.null(exactN)) exactN <- 0
                temp <- trappingtelemetry (
                           as.double(df0),
                           as.double(sigma),
                           as.double(z),
                           as.matrix(animals),
                           as.integer(noccasions),
                           as.integer(detectfn),
                           as.double(truncate^2),
                           as.integer(detectpar$binomN),
                           as.integer(exactN),
                           as.integer(maxperpoly)
                           )
            }
            if (temp$resultcode != 0) stop("result code", temp$resultcode, "from external simulation code")
            w <- array(dim=c(noccasions, nk, temp$n),
                dimnames = list(1:noccasions, polynames, rownames(animals)[temp$caught]))
            
            if (temp$n > 0) {
                w[,,] <- temp$value[1:prod(dim(w))]
            }
            w <- aperm(w, c(3,1,2))
            class(w) <- 'capthist'  ## moved here so we can use animalID(w) 2017-01-11
            traps(w)             <- traps
            
            if (temp$resultcode != 0) {
                if (temp$resultcode == 2)
                    stop ("more than ", maxperpoly, "  detections per animal",
                          " per polygon per occasion")
                else
                    stop ("call to ",simfunctionname, " failed")
            }
            detectedXY <- NULL
            if (temp$n > 0) {
                ## put XY coordinates in attribute
                nd <- sum(abs(w))
                detectedXY <- data.frame(matrix(ncol = 2,
                    temp$detectedXY[1:(2*nd)]))
                names(detectedXY) <- c('x','y')
                
            }
            if (detector(traps)[1] == 'telemetry') { 
                ID <- animalID(w)
                xyl <- split(detectedXY, ID)
                telemetryxy(w) <- xyl
            }
            else {
                xy(w) <- detectedXY
            }
        }
        ##-----------------------------------------------------------------------
        else if (detector(traps)[1] %in% c('signal','signalnoise')) {
            if (any(detectpar$binomN != 1))
                stop ("binomN != 1 not yet implemented for signal detectors")
            temp <- trappingsignal (
                as.double(beta0),
                as.double(beta1),
                as.double(sdS),
                as.double(cutval),
                as.double(muN),       # used only if signalnoise detector type
                as.double(sdN),       # used only if signalnoise detector type
                as.double(sdM),
                as.matrix(animals),
                as.matrix(traps),
                as.matrix(distmat2),
                as.matrix(usge),
                as.integer(detectfn)
            )
            # n = integer(1),
            # caught = integer(N),
            # signal = double(N*noccasions*k),
            # noise = double(N*noccasions*k),
            # value = integer(N*noccasions*k),    # detected/not detected
            # resultcode = integer(1)             # 0,1,2
            
            if (temp$resultcode != 0)
                stop ("call to 'trappingsignal' failed with resultcode ", temp$resultcode)
            w <- array(dim=c(noccasions, k, temp$n), dimnames =
                           list(1:noccasions,NULL,NULL))
            if (temp$n>0)  {
                w[,,] <- temp$value[1:(temp$n * noccasions * k)]
            }
            w <- aperm(w, c(3,1,2))
            class(w) <- 'capthist'
            traps(w) <- traps
        }
        else stop ('Unrecognised detector type')
        ##-----------------------------------------------------------------------

        if (!is.null(covariates(popn))) {
            if (all(detector(traps) %in% c('single','multi','proximity','count','capped'))) 
                covariates(w) <- covariates(popn)[row.names(w),, drop=F]
            else
                covariates(w) <- covariates(popn)[as.logical(temp$caught),, drop=F]
        }
        
        attr(w, 'cutval')    <- cutval
        attr(w, 'seed')      <- RNGstate      ## save random seed
        attr(w, 'detectpar') <- detectpar
        session(w)           <- '1'           ## dummy session values for now

        ## optionally save population, whether simulated or input  2014-04-27
        if (savepopn)
            attr(w, 'popn') <- popn

        if (detector(traps)[1] %in% c('signal','signalnoise')) {
            if (temp$n>0)  {
                signal(w) <- temp$signal[1:sum(w)]
                if (detector(traps)[1] %in% c('signalnoise'))
                    noise(w) <- temp$noise[1:sum(w)]
            }
        }

        if (renumber & (nrow(w)>0)) 
            rownames(w) <- 1:nrow(w)
        else {
            if (!all(detector(traps) %in% c('single','multi','proximity','count','capped')) && !HPXpoly) {
                rown <- rownames(popn)[temp$caught > 0]
                caught <- temp$caught[temp$caught>0]
                rownames(w) <- rown[order(caught)]
            }
        }

        w
    }
}
############################################################################################

sim.resight <- function (traps, popn = list(D = 5, buffer = 100, Ndist = 'poisson'), ...,
                         pID = 1, unmarked = TRUE, nonID = TRUE, unresolved = FALSE, 
                         unsighted = TRUE, pmark = 0.5, Nmark = NULL, markingmask = NULL) {

    ## markocc is vector length = noccasions; 0 = resighting, 1 = marking
    ## pmark is probability an individual is marked, across the entire population,
    ##     if all occasions are sighting occasions
    ############################################################################
    ## checks & setup

    if (is.null(markocc(traps))) {
        warning ("using default markocc 1 0 0 0 0")
        markocc <- c(1,0,0,0,0)
    }
    else
        markocc <- markocc(traps)

    markocc(traps) <- NULL  ## discard!
    S <- length(markocc)    ## number of sampling occasions
    K <- ndetector(traps)
    allsighting <- !any(markocc>0)
    unres <- markocc == -1
    proximityocc <- any(detector(traps) %in% 'proximity')

    dettype <- detectorcode(traps, MLonly = FALSE, noccasions = S)
    if (!all(dettype[markocc<1] %in% c(1,2,6,7)))
        stop ("only for sightings at binary or count proximity detectors", 
              " or polygon or transect searches")

    ## special case: Nmark animals distributed according to mask covariate
    distributedmarking <- FALSE
    if (!is.null(markingmask) & !is.null(Nmark)) {
        if (!is.null(covariates(markingmask)))
            distributedmarking <- 'marking' %in% names(covariates(markingmask))
    }
    
    ############################################################################
    ## make complete, identified capthist including all detections
    ## marking and sighting

    dots <- list(...)
    dots$traps <- traps

    #################################################################
    ## 2017-07-26
    seed <- dots$seed
    dots$seed <- NULL
    
    ## set random seed
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
    #################################################################
    
    dots$noccasions <- S    ## override noccasions
    dots$renumber <- FALSE  ## to match animalID

    if (!inherits(popn,'popn'))         ## generate popn if not provided
    {
        popn <- replacedefaults(list(D = 5, buffer = 100,
                                     Ndist = 'poisson'), popn)
        if (!(distributedmarking & allsighting))
            popn <- sim.popn (popn$D, core = traps, buffer = popn$buffer,
                              covariates = NULL, Ndist = popn$Ndist)
    }
    #---------------------------------------------------------------------------
    if (allsighting) {
        if (distributedmarking) {
           
            ## 2017-05-26
            ## Nunmark <- if (popn$Ndist == 'fixed') {
            Ndist <- attr(popn, 'Ndist')
            if (is.null(Ndist)) {
                warning ("popn does not have Ndist attribute; assuming Poisson")
                Ndist <- "poisson"
            }
            Nbuffer <- attr(popn, 'Nbuffer')
            D <- attr(popn, 'D')
            if (is.null(Nbuffer) & (is.null(D) | inherits(D, 'mask')))
                stop ("sim.resight option not available for missing Nbuffer,D or D as mask")
            
            Nunmark <- if (Ndist == 'fixed') {
                    if (!is.null(Nbuffer)) Nbuffer - Nmark
                else D * masksize(markingmask) - Nmark
            }
            else { 
                rpois(1, D * masksize(markingmask)) - Nmark
            }
            if (Nunmark<0) {
                warning ("Nunmark < 0; setting to 0")
                Nunmark <- 0
            }
    
            covariates(markingmask)$pi.marking <- covariates(markingmask)$marking / sum(covariates(markingmask)$marking)
            covariates(markingmask)$pi.unmarking <- 1 - covariates(markingmask)$marking
            popnM <- sim.popn(D = 'pi.marking', core = markingmask, Ndist = 'fixed', Nbuffer = Nmark, model2D='IHP')
            popnU <- sim.popn(D = 'pi.unmarking', core = markingmask, Ndist = 'fixed', Nbuffer = Nunmark, model2D='IHP')
        }
        else {
            gotten <- rep(TRUE, nrow(popn))
            
            ## circumscribe marked animals
            if (!is.null(markingmask)) {
                gotten <- pointsInPolygon(popn, markingmask)
                if (!is.null(covariates(markingmask)))
                    if ('marking' %in% names(covariates(markingmask)))
                        warning("uniform distribution over marking mask;",
                                " covariate 'marking' ignored")
            }
            
            if (!is.null(Nmark)) {
                if (Nmark > sum(gotten)) {
                    warning("Nmark exceeds population size, marking all")
                }
                else {
                    gottenN <- sample.int(sum(gotten), size = Nmark, replace = FALSE)
                    gotten[gotten] <- (1:sum(gotten)) %in% gottenN
                }
            }
            else {
                gotten <- gotten & (runif(length(gotten)) < pmark)
            }
            popnM <- subset(popn, gotten)     ## inside A0 and marked
            popnU <- subset(popn, !gotten)    ## outside A0 or not marked
        }
        dots$popn <- popnM
    }   ## end of allsighting block
    else {
        dots$popn <- popn
    }
    #---------------------------------------------------------------------------
    if (ms(popn) & length(popn)==2) {
        ## experimental combination of 'marking' and 'sighting' populations
        ## these have distinct x-y, but rows corresp to the same individuals
        lambda0 <- dots$detectpar$lambda0 
        usge <- usage(dots$traps)
        dots$popn <- popn[[1]]
        detector(dots$traps) <- detector(traps)[markocc==1]
        markocc(dots$traps) <- NULL
        dots$noccasions <- length(detector(dots$traps))
        usage(dots$traps) <- usage(dots$traps)[,markocc==1]
        dots$detectpar$lambda0 <- lambda0[markocc==1]
        CH1 <- do.call('sim.capthist', dots)
        
        dots$popn <- popn[[2]]
        detector(dots$traps) <- detector(traps)[markocc==0]
        dots$noccasions <- length(detector(dots$traps))
        usage(dots$traps) <- usge[,markocc==0]
        dots$detectpar$lambda0 <- lambda0[markocc==0]
        CH2 <- do.call('sim.capthist', dots)
        
        ## coerce detectors purely for the sake of join()
        detector(traps(CH1)) <- 'count'
        detector(traps(CH2)) <- 'count'
        capthist <- join(list(CH1, CH2))
        detector(traps(capthist)) <- detector(traps)
    }
    else {
        capthist <- do.call('sim.capthist', dots)
    }
    ############################################################################
    ## transform simulated capthist object into resighting data

    ## 2021-05-19 using ksn as safe universal order regardless of polygon/signal
    df <- data.frame(session  = rep(1,sum(abs(capthist))),
                     animalID = animalID(capthist, sortorder = 'ksn'), 
                     occasion = occasion(capthist, sortorder = 'ksn'),
                     trapID   = trap(capthist, sortorder = 'ksn'), 
                     stringsAsFactors = FALSE)
    if (all(detector(traps) %in% c('polygon', 'transect'))) {
        df <- cbind(df[,-4], xy(capthist))
        df$trapID <- xyinpoly(df[,4:5], traps(capthist))
    }
    df <- df[order(df$animalID, df$occasion, df$trapID),]

    if (nrow(df) == 0) {
        dfID <- data.frame (session = 1, animalID = "NONE", occasion = 1, trapID = 1)
        dfnonID <- dfID[-1,]
    }
    else {
        # data frame of detections on marking occasions
        dfmarking <- df[markocc[df$occasion]>0,,drop=FALSE]
        # when was each animal first marked?
        df$markingoccasion <- dfmarking$occasion[match(df$animalID, dfmarking$animalID)]
        df$markingoccasion[is.na(df$markingoccasion)] <- S+1    ## never marked
        if (allsighting) df$markingoccasion <- rep(0, nrow(df)) ## all marked
        # separate data frames for detections of marked and unmarked animals
        dfmarked <- df[df$occasion >= df$markingoccasion,, drop = FALSE]
        dfunmarked <- df[df$occasion < df$markingoccasion,, drop = FALSE]
        # which of the marked detections were identified?
        markingoccasions <- (1:S)[markocc>0]
        dfmarked$ID <- (dfmarked$occasion %in% markingoccasions) |
            (runif(nrow(dfmarked)) < pID)
        # separate data frames for detections of ID and nonID marked animals
        dfnonID <- dfmarked[!dfmarked$ID, , drop = FALSE]
        dfID    <- dfmarked[dfmarked$ID, , drop = FALSE]
    }

    ############################################################################
    ## Now build the new capthist object from dfID
    if (all(detector(traps) %in% c('polygon', 'transect'))) {
        CH <- make.capthist(dfID[,1:5, drop=FALSE], traps, fmt = 'XY', noccasions = S)
    }
    else {
        CH <- make.capthist(dfID[,1:4, drop=FALSE], traps, fmt = 'trapID', noccasions = S)
    }
    markocc(traps(CH)) <- markocc
    if (all(detector(traps(capthist)) %in% c('polygon', 'transect')))
        trapnames <- levels(polyID(traps(CH)))
    else
        trapnames <- row.names(traps(CH))
    
    ## 2017-09-01
    CH[,unres,] <- 0
    
    ############################################################################
    ## add all-zero histories for animals marked but never sighted (allsighting only)
    if (allsighting & unsighted) {
        covariates(CH) <- NULL  # beats problem in addzeroCH 2016-11-19
        nzero <- nrow(popnM) - nrow(CH)
        CH <- addzeroCH(CH, nzero)
    }
    ############################################################################
    ## add sightings of unmarked animals
    ## apply same sampling to unmarked fraction of population
    if (allsighting) {
        dots$popn <- popnU
        capthistU <- do.call('sim.capthist', dots)
    }
    
    ##------------------------------------------------------------------------------
    ## sightings of definitely unmarked animals
    
    if (unmarked) {
        countfn <- function(x) {
            x <- x * col(x)   ## x 0/1
            tabulate(x[x>0], nbins = K)
        }
        
        if (allsighting) {
            Tu <- as.matrix(table(
              factor(trap(capthistU, sortorder = 'snk'), levels = trapnames),
              factor(occasion(capthistU, sortorder = 'snk'), levels = 1:S)))
        }
        else {
            Tu <- as.matrix(table(factor(dfunmarked$trapID, levels = trapnames),
                                  factor(dfunmarked$occasion, levels = 1:S)))
        }
                
        ## 2017-03-16, presence only if detector == proximity
        Tu[, proximityocc] <- pmin (1, Tu[, proximityocc])
        Tu[,markocc>0] <- 0  ## just in case...
        Tu[,unres] <- 0
        Tu(CH) <- Tu
    }
    
    ##------------------------------------------------------------------------------
    ## sightings of marked animals that were not individuated
    
    if (nonID) {
        Tm <- as.matrix(table(factor(dfnonID$trapID, levels = trapnames),
              factor(dfnonID$occasion, levels = 1:S)))
        Tm[, proximityocc] <- pmin (1, Tm[, proximityocc])
        Tm[,unres] <- 0
        Tm(CH) <- Tm 
    }
    
    ##------------------------------------------------------------------------------
    ## sightings on occasions when mark status not recorded
    
    if (unresolved) {
        if (any(unres)) {
            Tn <- as.matrix(table(factor(df$trapID, levels = trapnames),
                                  factor(df$occasion, levels = 1:S)))
            + as.matrix(table(factor(dfunmarked$trapID, levels = trapnames),  
                              factor(dfunmarked$occasion, levels = 1:S)))
            Tn[,!unres] <- 0
            Tn(CH) <- Tn
            ## optionally drop all-zero histories
            if (allsighting & !unsighted) {
                OK <- apply(CH,1,sum)>0
                CH <- subset(CH,OK)
            }
        }    
    }
    
    ##------------------------------------------------------------------------------
    
    ## if savepopn = TRUE...
    if (!is.null(attr(capthist, 'popn'))) {
        if (allsighting) {
            popn <- rbind(popnM, popnU)
            covariates(popn) <- data.frame(marked = rep(c(TRUE, FALSE), c(nrow(popnM), nrow(popnU))))
            attr(CH, 'popn') <- popn
        }
        else {
            ## 2015-12-18; this could be moved to sim.capthist
            covariates(popn) <- data.frame(marked = rownames(popn) %in% dfmarked$animalID)   
            attr(CH, 'popn') <- popn   ## NULL unless 
        }
    }
    attr(CH, 'seed') <- RNGstate      ## save random seed
    CH
}
############################################################################################

