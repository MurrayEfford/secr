############################################################################################
## package 'secr'
## autoini.R
############################################################################################

autoini <- function (capthist, mask, detectfn = 0, thin = 0.2, tol = 0.001,
                     binomN = 1, adjustg0 = TRUE, adjustsigma = 1.2, ignoreusage = FALSE, 
                     ncores = NULL)

# obtain approximate fit of HN SECR model
# for use as starting values in MLE
# uses external C code
# binomN is scalar
    
{
    naivedcall <- function (sigma)
    {
      db - naivedcpp(
        sigma,
        wt,
        as.matrix(trps),
        as.matrix(mask),
        detectfn)
    }
    
    naivecap3 <- function (lambda0, sigma, cap)
    {
        cap - naivecap3cpp(dettype[1], 
                           lambda0, 
                           sigma, 
                           as.matrix(usge),
                           as.matrix(trps), 
                           as.matrix(mask),
                           14)
    }
    
    naiveesa <- function (g0, sigma)
    {
      nc <- 1
      m <- nrow(mask)
      realparval0 <- matrix(c(g0,sigma), nrow=1, ncol=2)
      PIA0 <- array(1, dim=c(1,nc,s,k,1))
      # allow for binary use/non-use of detectors
      if (!is.null(usage(trps))) PIA0[t(usage(trps)==0)] <- -1
      
      distmat2 <- secr_getdistmat2(trps, mask, NULL)
      gkhk <- makegkPointcpp (
          as.integer(detectfn), 
          as.integer(grain), 
          as.integer(ncores),
          as.matrix(realparval0), 
          as.matrix(distmat2), 
          as.double(0))
      
      if (any(dettype==0)) CH0 <- secr_nullCH (c(nc,s), FALSE)
      else CH0 <- secr_nullCH (c(nc,s,k), FALSE)
      binomNcode <- secr_recodebinomN(dettype, binomN, 0)
      pmixn <- matrix(1, nrow=1, ncol=nc)
      pdot <- secr_integralprw1 (
          cc0 = nrow(realparval0), 
          haztemp = secr_gethazard(m, binomNcode, nrow(realparval0), gkhk$hk, PIA0, usge), 
          gkhk = gkhk, 
          pi.density = matrix(1/m, nrow=m, ncol=1), 
          PIA0 = PIA0, 
          CH0 = CH0, 
          binomNcode = binomNcode, 
          MRdata = list(markocc = rep(1,s), firstocc = -1),  ## all marking
          grp = rep(1,nc), 
          usge = usge, 
          pmixn = pmixn,
          pID = matrix(1, nrow=s, ncol=1),
          grain = as.integer(grain),
          ncores = ncores)
      pdot * masksize(mask)
    }
    ##########################
    ## main line
    computeD <- TRUE
    if (length(tol)==1) tol <- rep(tol,2)

    if (nrow(capthist)<5)
        stop ("too few values for autoini")  ## message changed 2015-01-06

    if (is.character(detectfn))
        detectfn <- secr_detectionfunctionnumber(detectfn)

    if (! (detectfn %in% c(0)))
        stop ("only halfnormal detection function implemented in 'autoini'")

    ncores <- setNumThreads(ncores)
    grain <- if (!is.null(ncores) && (ncores==1)) 0 else 1
    obsRPSV <- RPSV(capthist, CC = TRUE)

    # drop telemetry component
    ## drop zero histories
    ## drop occasions with dettype==13
    ## drop arbitrary telemetry detector
    # ncapt <- apply(capthist,1,sum)
    # 
    # if (any(ncapt==0) || any (dettype==13)) {
    #     nontelemetry <- dettype != 13
    #     capthist <- subset(capthist, ncapt>0, occasions = nontelemetry)   ## drop zero histories
    #     binomN <- binomN[nontelemetry]
    #     dettype <- dettype[nontelemetry]
    # }
    
    trps <- traps(capthist)
    usge <- usage(trps)
    if (is.null(usge) || ignoreusage) {
        ## assuming k = nk i.e. not polygon or transect detector
        usge <- matrix(1, nrow = nrow(trps), ncol = ncol(capthist))
    }
    dettype <- secr_detectorcode(trps, noccasions = ncol(capthist))

    n       <- nrow(capthist)    # number of individuals
    s       <- ncol(capthist)    # number of occasions
    k       <- nrow(trps)
    markocc <- markocc(traps(capthist))
    if (is.null(markocc)) 
        markocc <- rep(1,s)     ## assume all occasions were simple marking occasions
    allsighting <- !any(markocc>0)

    if (!all(dettype %in% c(-1:5,8,13)))
        list(D=NA, g0=NA, sigma=NA)
    else {
        
        ## treat capped as proximity here 2022-04-13
        dettype[dettype==8] <- 1  
        
        ## wt is the number of opportunities for capture given binary usage
        wt <- apply(usge>0, 1, sum)

        # optionally thin mask
        if ((nrow(mask)>100) & (thin>0) & (thin < 1))
            mask <- mask[runif(nrow(mask)) < thin,]
        else
            thin <- 1.0
        m        <- nrow(mask)                   # number of points in mask
        cpa     <- sum(abs(capthist))/n      # captures per animal
        
        if (is.na(obsRPSV) | (obsRPSV<1e-10)) {    ## try db
            db <- dbar(capthist)
            if (!is.null(attr(trps,'spacing'))) {
                if (is.na(db)) {
                    warning ("could not calculate 'dbar'; using detector spacing")
                    db <- attr(trps, 'spacing')
                }
                if (db < (attr(trps, 'spacing')/100)) {
                    warning ("'dbar' close to zero; using detector spacing instead")
                    db <- attr(trps, 'spacing')
                }
            }
            if (is.na(db) | is.nan(db) | (db<1e10) )
                return(list(D=NA,g0=NA,sigma=NA))
            else
                tempsigma <- uniroot (naivedcall, lower = db/10, upper = db*10, tol=tol[2])$root
        }
        else {
            tempsigma <- obsRPSV * adjustsigma  
        }
        if (is.null(usage(trps))) wt <- rep(s,k)
        low <- naivecap3(0.00001, sigma=tempsigma, cap=cpa)
        upp <- naivecap3(100, sigma=tempsigma, cap=cpa)
        badinput <- FALSE
        if (is.na(low) | is.na(upp)) badinput <- TRUE
        else if (sign(low) == sign(upp)) badinput <- TRUE
        
        if (badinput) {
            # not sure what conditions cause this 28/4/2008
            # observed number cap more than expected when g0=1 28/8/2010
            # maybe better in future to set g0 = 0.9
            warning ("'autoini' failed to find g0; setting initial g0 = 0.1")
            tempg0 <- 0.1
        }
        else {
            templambda0 <- uniroot (naivecap3, lower=0.00001, upper=0.99999,
                 f.lower = low, f.upper = upp, tol=tol[1],
                 sigma=tempsigma, cap=cpa)$root
            tempg0 <- 1 - exp(-templambda0)
        }
        
        if (computeD) {
            if (allsighting)  ## includes all-zero rows
                tempD <- n / masksize(mask)
            else {
                esa <- naiveesa (tempg0, tempsigma)
                tempD <-  n / esa * thin
            }
        }
        else tempD <- NA

        ## 2012-12-18,24 adjust for large effort and/or binomN
        ## when is this needed? 2019-09-03
        if (adjustg0) {
            if (binomN == 0)
                adjusted.g0 <- tempg0 / usge[usge>0]
            else
                adjusted.g0 <- 1 - (1 - tempg0)^ ( 1 / (usge[usge>0] * binomN) )
            tempg0 <- mean(adjusted.g0)
        }

        list(D = tempD, g0 = tempg0, sigma = tempsigma)
    }
}
##################################################################################
