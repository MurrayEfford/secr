###############################################################################
## package 'secr'
## preparedata.R
###############################################################################

#--------------------------------------------------------------------------------
getsignal <- function (dettype, capthist, tx) {
  if (all(dettype %in% c(5,12))) {    # signal strength, signalnoise
    signal <- signalmatrix(capthist)
    noise <- signalmatrix(capthist, noise = TRUE)
    signal <- switch( tx,
                      log = log(signal),
                      logit = logit(signal),
                      identity = signal
    )
    noise <- switch( tx,
                      log = log(noise),
                      logit = logit(noise),
                      identity = noise
    )
    signal[is.na(signal)] <- -1   ## for C++ code
    noise[is.na(noise)] <- -1     ## for C++ code
    list(signal = signal, noise = noise)
  }
  else {
    list(signal=-1, noise=-1)
  }
}
#--------------------------------------------------------------------------------

secr_recodebinomN <- function (dettype, binomN, telemcode) {
  binomN <- secr_expandbinomN(binomN, dettype)
  detectr <- .localstuff$validdetectors[dettype+2]
  detectr[(detectr %in% c('count','polygon', 'transect')) & (binomN == 0)] <- "poissoncount"
  detectr[(detectr %in% c('count','polygon', 'transect')) & (binomN > 0)]  <- "binomialcount"
  newbinomN <- function (det,N) {
      recoded <- switch(det, 
           single = -2, multi = -2, polygonX = -2, transectX = -2,
           proximity = -1, signal = -1, capped = -1, 
           poissoncount = 0, binomialcount = N, telemetry = -3, -9)
    ## dummy value -9 indicates no action
    if (recoded == -3 && telemcode == 0) recoded <- -7
    recoded
  }
  out <- mapply(newbinomN, detectr, binomN)
  if (any(out < -9))
    stop("secr not ready for detector type")
  out
}
#--------------------------------------------------------------------------------

secr_nullCH <- function (dimCH, individual) {
    if (is.null(individual)) {
        individual <- TRUE   ## 2020-05-16 for backward compatibility
    }
    if (!individual) {
        dimCH[1] <- 1
    }
    array(0, dim = dimCH)
}
#--------------------------------------------------------------------------------

# compress ch (n x k) by animal: list of detectors with non-zero counts, terminated by -1
# (n x k x 2)
nk2 <- function(ch) {
    one <- function(i) {
        wh <- which(ch[i,1,]>0)
        out <- matrix(-1, nrow = nk, ncol = 2)
        if (length(wh)>0) {
            out[1:length(wh),] <- cbind(ch[i,1,wh], wh-1)
        }
        out
    }
    nk <- dim(ch)[3]
    ## 2022-01-04 catch zero rows
    if (nrow(ch) < 1) {
      array(dim=c(0,nk,2))
    }
    else {
      ch2 <- sapply(1:nrow(ch), one)
      ch3 <- array(ch2, dim = c(nk,2,nrow(ch)))
      aperm(ch3, c(3,1,2))
    }
}
#--------------------------------------------------------------------------------

## 2-D n x s if exclusive detector
## 2-D s x k if capped detector (deferred)
## 3-D n x k x 2 (which detector, count of positive records)
compressCH <- function (CH, binomN, fastproximity) {  
    if (all(binomN == -2)) {
        lost <- apply(CH, 1:2, min)<0
        CH <- abs(CH)
        CH <- apply(CH, 1:2, which.max) *(apply(CH, 1:2, max)>0)
        CH[lost] <- -CH[lost]
    }
    # else if (all(binomN == -3)) {
    #   lost <- apply(CH, 2:3, min)<0
    #   CH <- abs(CH)
    #   CH <- apply(CH, 2:3, which.max) *(apply(CH, 2:3, max)>0)
    #   CH[lost] <- -CH[lost]
    # }
    else if (fastproximity) {
        CH <- nk2(CH) 
    }
    CH
}
############################################################################################

decompressCH <- function (CH, fastproximity) {  
    if (fastproximity) {
        out <- array(0, dim=c(nrow(CH), 1, ncol(CH)))
        for (i in 1:nrow(CH)) {
            n <- which(CH[i,,1]>0)
            k <- CH[i,n,1]
            count <- CH[i,n,2]
            out[i,1,k] <- count
        }
        return(out)
    }
    else {
        return (CH)
    }
}
############################################################################################

## settings for mark-resight
markresightdata <- function (capthist, mask, fixed, chat, control, knownmarks) {
    getsight <- function(T) {
        Tval <- attr(capthist,T)
        tmp <- if ((control[[T]]=='ignore') | is.null(Tval))
            NULL
        else
            if (control[[T]]=='sum')
                sum(Tval)
        else
            if (control[[T]]=='bydetector') {
                if (is.matrix(Tval)) {
                    apply(Tval, 1, sum)
                }
                else {
                    if (length(Tval) == secr_ndetector(traps(capthist)))
                        Tval
                    else
                        stop ("bydetector expects", T, "as a matrix or length-K vector,",
                              "where K is the number of detectors")
                }
            }
        else
            Tval
        ## second element is number of values; '1' indicates summed counts
        if (is.null(tmp)) NULL else tmp
    }
    markocc <- markocc(traps(capthist))
    s <- ncol(capthist)
    if (is.null(markocc)) {
        markocc <- rep(1, s)
        Tu <- Tm <- Tn <- NULL
        allsighting <- FALSE
        anysighting <- FALSE
        firstocc <- rep(-1,nrow(capthist))
    }
    else {
        m <- nrow(mask)
        defaultcontrol <- list(Tu='as.is', Tm='as.is', Tn='ignore')
        # possible control values
        #   ignore
        #   as.is
        #   bydetector
        #   sum
        control <- secr_replacedefaults(defaultcontrol, control)
        allsighting <- !any(markocc>0)
        anysighting <- any(markocc<1)
        
        ## if (CL) control$Tu <- 'ignore'
        if (!any(markocc==0)) control$Tm <- 'ignore'
        
        
        ## risk of double counting:
        ## consider Tn only when markocc[s] = -1
        ## there should be no Tu on those occasions
        if (any(markocc<0)) control$Tn <- 'as.is'
        
        if (!is.null(fixed$pID)) {
            if (fixed$pID == 1) control$Tm <- 'ignore'
        }
        if(is.null(fixed$pID) & control$Tm == 'ignore')
            warning("Set fixed = list(pID=1) if no sightings of unidentified marked animals Tm")
        
        Tu <- getsight('Tu')
        Tm <- getsight('Tm')
        Tn <- getsight('Tn')
        
        if (allsighting) {
            ## assume all to be pre-marked & available for detection
            firstocc <- rep(-1,nrow(capthist))
        }
        else {
            ch2 <- apply(abs(capthist),1:2,sum)
            firstocc <- apply(ch2>0,1,match, x=TRUE)-1
            firstocc[is.na(firstocc)] <- s
        }
        
        ## special case: unmarked or presence/absence detector
        detect <- detector(traps(capthist))
        if (any(detect %in% c('unmarked','markocc'))) {
            # no action needed?
        }
    }
    if (!is.null(chat)) {
        # if (is.matrix(chat))
        #     chat <- chat[sessnum,]
        # else {
        chat <- unlist(chat)
        if (length(chat)==1) {
            chat <- c(chat,1,chat)
            warning("assuming chat 1.0 for Tm")
        }
        if (any(chat<1)) {
            warning("setting chat < 1.0 to 1.0")
        }
        chat <- pmax(chat, 1)
        ## }
    }
    else
        chat <- c(1,1,1)
    
    
    pi.mask <- -1      ## signals pimask not used
    sightmodel <- 0
    if (allsighting) {
        ## pi.mask is Pr(marked animal is from pixel m)
        ## i.e. pdf(x) * area
        pi.mask <- rep(1/nrow(mask), nrow(mask))
        if (!is.null(maskcov <- covariates(mask))) {
            if ('marking' %in% names (maskcov)) {
                if (any(is.na(maskcov$marking)) | any (maskcov$marking<0))
                    stop ("invalid marking covariate in mask")
                pi.mask <- maskcov$marking / sum (maskcov$marking)
            }
        }
        if (knownmarks)
            sightmodel <- 5
        else 
            sightmodel <- 6
    }
    
    list(markocc = markocc, Tu = Tu, Tm = Tm, Tn = Tn,
         anysighting = anysighting, allsighting = allsighting,
         chat = chat, pi.mask = pi.mask, firstocc = firstocc, 
         sightmodel = sightmodel)
}

##############################################################################

secr_prepareSessionData <- function (capthist, mask, maskusage, 
    design, design0, detectfn, groups, fixed, hcov, details, 
    aslist = TRUE, sessnum = 1) {
    ## aslist used internally to determine whether single-session data are wrapped in a list
    if (ms(capthist)) {
        if (!ms(mask)) stop ("expect session-specific mask in secr_prepareSessionData")
        if (is.null(maskusage))
            maskusage <- vector('list', length(capthist))
        mapply(
            secr_prepareSessionData, 
            capthist = capthist, 
            mask = mask, 
            maskusage = maskusage,
            sessnum = 1:length(capthist),
            MoreArgs = list(design, design0, detectfn, groups, fixed, 
                hcov, details, FALSE), 
            SIMPLIFY = FALSE)
    }
    else {
        nc   <- nrow(capthist)
        s    <- ncol(capthist)
        m    <- nrow(mask)
        traps   <- traps(capthist)
        dettype <- secr_detectorcode(traps, MLonly = TRUE, noccasions = s)
        binomNcode <- secr_recodebinomN(dettype, details$binomN, secr_telemcode(traps))
        ## k-1 because we have zero-terminated these vectors
        k <- secr_getk(traps)
        K <- if (length(k)>1) length(k)-1 else k
        cumk <- cumsum(c(0,k))[1:length(k)]
        
        ## mark-resight
        MRdata <- markresightdata(capthist, mask, fixed,
            details$chat, details$markresight, details$knownmarks)

        ## knownclass for hcov mixture models
        knownclass <- secr_getknownclass(capthist, details$nmix, hcov)
        
        ## get static distance matrix
        distmat2 <- secr_getdistmat2(traps, mask, details$userdist)

        n.distrib <- switch (tolower(details$distribution), poisson=0, binomial=1, 0)
        signal <- getsignal (dettype, capthist, details$tx)
        xy <- secr_getxy (dettype, capthist)
        usge <- usage(traps)
        if (is.null(usge) || details$ignoreusage) {
            usge <- matrix(1, nrow = K, ncol = s)
        }
        if (is.null(maskusage)) {
            maskusage <- secr_maskboolean(capthist, mask, details$maxdistance)
        }
        else {
            if (!is.matrix(maskusage) || nrow(maskusage) != nrow(capthist) || ncol(maskusage) != nrow(mask))
                stop ('specified maskusage should be n x m matrix of logical values')
            maskusage[] <- as.logical(maskusage)
        }
        
        if (!is.null(details$externalpdot)) {
            if (!(details$externalpdot %in% names(covariates(mask)))) 
                stop ("externalpdot '", details$externalpdot, "' not found in mask covariates")
            externalpdot <- covariates(mask)[, details$externalpdot]
            message("using external pdot")
        }
        else {
            externalpdot <- NULL
        }
        
        if (!is.null(details$externalqx)) {
            if (!(details$externalqx %in% names(covariates(mask)))) 
                stop ("externalqx '", details$externalqx, "' not found in mask covariates")
            externalqx <- covariates(mask)[, details$externalqx]
            if (!is.null(externalpdot)) stop ("specify only one of externalqx, externalpdot")
            message("using external q(x)")
        }
        else {
            externalqx <- NULL
        }
        
        ## Groups
        grp  <- secr_group.factor (capthist, groups)
        if (any(is.na(grp))) {
            stop("group is missing for at least one animal")
        }
        ngroup <- max(1,length(secr_group.levels(capthist, groups)))
        CH <- compressCH(capthist, binomNcode, details$fastproximity) 
        
        # 2023-06-09 tentatively remove hcov from condition
        # this leaves some uncertainty: 
        # when is full CH0 (1 row per animal) really needed?
        # why is this an issue for polygonhistoriescpp and not simplehistoriescpp?
        
        CH0 <- secr_nullCH(dim(CH), packageVersion('secr')<'4.0.0' || design0$individual || ngroup>1)   ## all-zero CH
        
        #####################################################################
        ## unclear whether this is correct wrt groups
        if (all(detector(traps) %in% .localstuff$simpledetectors)) {
            logmult <- logmultinom(capthist, secr_group.factor(capthist, groups))
        }
        else {
            logmult <- 0
        }
        #####################################################################
        
        data <- list(
            sessnum = sessnum,    # added 2021-06-22
            CH = CH,
            CH0 = CH0,
            nc = nc,
            s = s,
            k = k,
            K = K,
            cumk = cumk,
            m = m,
            traps = traps,
            dettype = dettype,
            binomNcode = binomNcode,
            usge = usge,
            mask = mask,
            externalpdot = externalpdot,
            externalqx = externalqx,
            distmat2 = distmat2,
            knownclass = knownclass,
            n.distrib = n.distrib,
            MRdata = MRdata,
            signal = signal,
            xy = xy,
            grp = grp,
            maskusage = maskusage,
            logmult = logmult
        )
    if (aslist) list(data=data)
    else data
    #####################################################################
  }
}
