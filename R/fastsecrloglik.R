###############################################################################
## package 'secr'
## fastsecrloglik.R
## likelihood evaluation functions
## 2020-10-11 knownclass bug when not all classes present in session fixed
###############################################################################

#--------------------------------------------------------------------------------
allhistfast <- function (realparval, gkhk, pi.density, PIA, 
                         nk2ch, usge, pmixn, maskusage,
                         grain, ncores, binomN, indiv) {
    nc <- dim(nk2ch)[1] # dim(PIA)[2]
    ## 2022-01-04
    if (nc<1) return(1)
    nmix <- dim(PIA)[5]
    m <- length(pi.density)
    sump <- numeric(nc)
    
    for (x in 1:nmix) {
      temp <- fasthistoriescpp(
        as.integer(m),
        as.integer(nc),
        as.integer(nrow(realparval)),
        as.integer(grain),
        as.integer(ncores),
        as.integer(binomN),
        as.logical(indiv),
        matrix(nk2ch[,,1], nrow=nc),
        matrix(nk2ch[,,2], nrow=nc),
        as.double (gkhk$gk),  ## precomputed probability 
        as.double (gkhk$hk),  ## precomputed hazard 
        as.double (pi.density),
        as.integer(PIA[1,,,,x]),
        as.integer(usge),
        as.matrix (maskusage))
      sump <- sump + pmixn[x,] * temp
    }
    sump
}
#--------------------------------------------------------------------------------

secr_integralprw1fast <- function (realparval0, gkhk, pi.density, PIA0, 
                              nk2ch0, usge, pmixn, grain, ncores, binomN, indiv) {
    nc <- dim(PIA0)[2]
    nr <- nrow(nk2ch0)
    nmix <- dim(PIA0)[5]
    m <- length(pi.density)
    sump <- numeric(nc)
    for (x in 1:nmix) {
        temp <- fasthistoriescpp(
            as.integer(m),
            as.integer(nr),    ## 1 
            as.integer(nrow(realparval0)),
          as.integer(grain),
          as.integer(ncores),
          as.integer(binomN),
            as.logical(indiv),
            matrix(nk2ch0[,,1], nrow = nr),
            matrix(nk2ch0[,,2], nrow = nr),
            as.double (gkhk$gk),        ## precomputed probability 
            as.double (gkhk$hk),        ## precomputed hazard
            as.double (pi.density),
            as.integer(PIA0[1,1:nr,,,x]),
            as.integer(usge),
            as.matrix (matrix(TRUE, nrow = nr, ncol = m)))
        sump <- sump + pmixn[x,1:nc] * (1 - temp)
    }
    sump
}

#######################################################################################
fastsecrloglikfn <- function (
    beta, 
    parindx, 
    link, 
    fixed, 
    designD, 
    designNE, 
    design, 
    design0, 
    CL, 
    detectfn,
    learnedresponse,
    sessionlevels,
    data,
    details,
    dig = 3, betaw = 10, neglik = TRUE)

# Return the negative log likelihood for spatial capture-recapture model

# Transformed parameter values (density, g0, sigma, z etc.) are passed in the vector 'beta'
# 'detectfn' is integer code for detection function
#    0 = halfnormal, 1 = hazard, 2 = exponential etc.
# 'CL' is logical for conditional (CL=T) vs full (CL=F) likelihood
# details$trace=T sends a one-line report to the screen

{
    #--------------------------------------------------------------------------------
    sessionLL <- function (data, sessnum, like = 0) {
        ## log likelihood for one session
        ## in multi-session case must get session-specific data from lists
        ###################################################################
        #---------------------------------------------------
      PIA <- design$PIA[sessnum, 1:data$nc, 1:data$s, 1:data$K, ,drop=FALSE]
        ## unmodelled beta parameters, if needed
        miscparm <- secr_getmiscparm(details$miscparm, detectfn, beta, parindx, details$cutval)
        D.modelled <- (!CL || details$relativeD) && is.null(fixed$D)
        
        density <- secr_getmaskpar(D.modelled, D, data$m, sessnum, details$unmash, 
                              attr(data$capthist, 'n.mash'))
        if (!D.modelled) {
            pi.density <- rep(1/data$m, data$m)  
        }
        else {
            if (!is.null(details$externalqx)) {
                density <- density * data$externalqx
            }
            pi.density <- density / sum(density)
        }
        #---------------------------------------------------
        ## allow for scaling of detection
        Dtemp <- if (D.modelled) mean(D[,1,sessnum]) else NA
        Xrealparval <- secr_reparameterize (realparval, detectfn, details,
                                       data$mask, data$traps, Dtemp, 1) # 1 was s! 2024-07-30
        ## check valid parameter values
        if (!all(is.finite(Xrealparval))) {
            cat ('beta vector :', beta, '\n')
            warning ("extreme 'beta' in 'fastsecrloglikfn' ",
                     "(try smaller stepmax in nlm Newton-Raphson?)")
            return (1e10)
        }
        ## DOES NOT ALLOW FOR GROUP VARIATION IN DENSITY
        ## more thoughts 2015-05-05
        ## could generalize by
        ## -- making Dtemp a vector of length equal rows in realparval
        ## -- matching either
        ##      first group (as before)
        ##      sum of all groups
        ##      own group [PROBLEM: locating group of each realparval row]
        ## in all cases density is the mean over mask points
        
        ## CHECK use of Dtemp in regionN.R, sim.secr.R
        ## PERHAPS for consistency make a function to construct Dtemp vector
        ## given mask, model, group matching rule (first, sum, own)
        
        ##--------------------------------------------------------------

        #####################################################################
        pmixn <- secr_getpmix (data$knownclass, PIA, Xrealparval)  ## membership prob by animal
        if (!is.null(details$userdist)) {    # changed from is.function() 2024-02-15
            distmat2 <- secr_getuserdist(data$traps, data$mask, details$userdist, sessnum, 
                                    NElist, density[,1], miscparm = miscparm, detectfn = detectfn)
        }
        else {
            distmat2 <- data$distmat2
        }
        
        ## precompute gk, hk for point detectors
        if (data$dettype[1] %in% c(0,1,2,5,8)) {
            gkhk <- makegkPointcpp (
                as.integer(detectfn),
                as.integer(details$grain),
                as.integer(details$ncores),
                as.matrix(Xrealparval),
                as.matrix(distmat2),
                as.double(miscparm))
            if (details$anycapped) {   ## capped adjustment
              gkhk <- cappedgkhkcpp (
                as.integer(nrow(Xrealparval)),
                as.integer(nrow(data$traps)),
                as.double(attr(data$mask, "area")),
                as.double(density[,1]),
                as.double(gkhk$gk), as.double(gkhk$hk))
            }
        }
        
        prw <- allhistfast (Xrealparval, gkhk, pi.density, PIA, 
          data$CH, data$usge, pmixn, data$maskusage, 
          details$grain, details$ncores, details$binomN, design$individual)
        if (!is.null(details$externalpdot)) {
            pdot <- rep(sum(data$externalpdot * pi.density), data$nc)
        }
        else {
            pdot <- secr_integralprw1fast (Xrealparval, gkhk, pi.density, PIA, 
                  data$CH0, data$usge, pmixn, details$grain, details$ncores, 
                  details$binomN, design$individual)
        }
        if (details$debug>2) browser()
        
        comp <- matrix(0, nrow = 6, ncol = 1)
        
        #----------------------------------------------------------------------
        
        comp[1,1] <- if (any(is.na(prw)) || any(prw<=0)) NA else sum(log(prw))
        
        #----------------------------------------------------------------------
        
        if (!data$nc<=0) {
            comp[2,1] <- if (any(is.na(pdot)) || any(pdot<=0)) NA else -sum(log(pdot)) 
        }
        
        #----------------------------------------------------------------------
        
        if (!CL) {
            N <- sum(density[,1]) * secr_getcellsize(data$mask)
            ## 2022-01-05 catch nc = 0
            meanpdot <- if (data$nc == 0) pdot else data$nc / sum(1/pdot)
            ## 2023-09-22
            if (data$n.distrib == 1 && .localstuff$iter == 0 && data$nc>N) {
                warning("distribution = 'binomial' ",
                        "but number detected n (", data$nc, 
                        ") exceeds initial value of N (", round(N,1), ")")
            }
            comp[3,1] <- switch (data$n.distrib+1,
                                 dpois(data$nc, N * meanpdot, log = TRUE),
                                 secr_lnbinomial (data$nc, N, meanpdot),
                                 NA)
        }
        #----------------------------------------------------------------------
        
        ## adjustment for mixture probabilities when class known
        known <- sum(data$knownclass>1)
        if (details$nmix>1 && known>0) {
            nb <- details$nmix + 1
            nm <- tabulate(data$knownclass, nbins = nb)
            pmix <- attr(pmixn, 'pmix')
            
            # for (x in 1:details$nmix) {
            #     # need group-specific pmix
            #     comp[4,1] <- comp[4,1] + nm[x+1] * log(pmix[x]) 
            # }
            
            ## 2022-01-16 bug fix
            # firstx <- match ((1:details$nmix)+1, data$knownclass)
            # tempsum <- sum(pdot[firstx] * pmix)
            # comp[4,1] <- sum(nm[-1] * log(pdot[firstx] * pmix / tempsum))
            ##

            ## 2022-10-25 bug fix
            firstx <- match ((1:details$nmix)+1, data$knownclass)
            pdpmix <- pdot[firstx] * pmix
            pdpmix <- pdpmix[!is.na(pdpmix)]
            comp[4,1] <- sum(nm[-1] * log(pdpmix / sum(pdpmix)))
            ##
            
        }
        #----------------------------------------------------------------------
        
        if (details$debug>=1) {
            comp <- apply(comp,1,sum)
            cat(comp[1], comp[2], comp[3], comp[4], comp[5], comp[6], data$logmult, '\n')
        }
        sum(comp) + data$logmult
        
    } ## end sessionLL
    
    ###############################################################################################
    ## Main line of fastsecrloglikfn
    nsession <- length(sessionlevels)

    #--------------------------------------------------------------------
    # Fixed beta
    pbeta <- beta   # save varying beta, for trace etc.
    beta <- secr_fullbeta(beta, details$fixedbeta)

    #--------------------------------------------------------------------
    # Detection parameters
    detparindx <- parindx[!(names(parindx) %in% .localstuff$spatialparameters)]
    detlink <- link[!(names(link) %in% .localstuff$spatialparameters)]
    realparval  <- secr_makerealparameters (design, beta, detparindx, detlink, fixed)
    
    #--------------------------------------------------------------------
    sessmask <- lapply(data, '[[', 'mask')
    grplevels <- unique(unlist(lapply(data, function(x) levels(x$grp))))
    #--------------------------------------------------------------------
    # Density
    D.modelled <- (!CL || details$relativeD) && is.null(fixed$D)
    if (D.modelled) {
        D <- secr_getD (designD, beta, sessmask, parindx, link, fixed,
                   grplevels, sessionlevels, parameter = 'D')
        if (!is.na(sumD <- sum(D)))
            if (sumD <= 0)
                warning ("invalid density <= 0")
    }
    #--------------------------------------------------------------------
    # Non-Euclidean distance parameters
    NElist <- mapply(secr_getD, designNE, parameter = names(designNE),
                 MoreArgs = list(beta, sessmask, parindx, link, fixed,
                                 grplevels, sessionlevels), SIMPLIFY = FALSE)
    #--------------------------------------------------------------------
    # typical likelihood evaluation
    if (!is.null(details$saveprogress) && details$saveprogress>0 &&
        .localstuff$iter == 0) {
        secr_saveprogress(pbeta, NA, details$progressfilename)
    }
    
    loglik <- sum(mapply (sessionLL, data, 1:nsession))
    .localstuff$iter <- .localstuff$iter + 1   ## moved outside loop 2011-09-28
    if (details$trace) {
        cat(format(.localstuff$iter, width=4),
            formatC(round(loglik,dig), format='f', digits=dig, width=10),
            formatC(pbeta, format='f', digits=dig+1, width=betaw),
            '\n')
        flush.console()
    }
    #--------------------------------------------------------------------
    if (!is.null(details$saveprogress) && details$saveprogress>0 &&
        .localstuff$iter %% details$saveprogress == 0) {
        secr_saveprogress(pbeta, loglik, details$progressfilename)
    }
    
    #--------------------------------------------------------------------
    
    loglik <- ifelse(is.finite(loglik), loglik, -1e10)
    ifelse (neglik, -loglik, loglik)
}  ## end of fastsecrloglikfn
############################################################################################

