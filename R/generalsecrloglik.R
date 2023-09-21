###############################################################################
## package 'secr'
## secrloglik2.R
## likelihood evaluation functions
## 2019-10-12 moved helper fn to separate file
## 2019-12-04 integralprw1 modified to allow for individual covariate
## 2020-04-24 learnedresponse bug with multicatch traps fixed
## 2020-10-11 knownclass bug when not all classes present in session fixed
## 2021-08-07 ncores passed to C++ for parallelFor
###############################################################################

# dettype
# 
# single      = -1,
# multi       = 0,
# proximity   = 1,
# count       = 2,
# polygonX    = 3,
# transectX   = 4,
# signal      = 5,
# polygon     = 6,
# transect    = 7,
# capped      = 8,
# unmarked    = 10,
# presence    = 11,
# signalnoise = 12,
# telemetry   = 13,
# index       = 14

#--------------------------------------------------------------------------------
allhistsimple <- function (cc, haztemp, gkhk, pi.density, PIA, 
                           CH, binomNcode, MRdata, grp, usge, pmixn, pID, maskusage,
                           telemhr = 0, telemstart = 0, telemscale = 0,
                           grain, ncores, R = FALSE, debug = FALSE) {
  nc <- nrow(CH)
  ## 2022-01-04
  if (nc<1) return(1)
  k <- nrow(usge)
  m <- nrow(pi.density)
  nmix <- nrow(pmixn)
  ngroup <- length(unique(grp))
  sump <- numeric(nc)
  
  for (x in 1:nmix) {
      hx <- if (any(binomNcode==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## lookup sum_k (hazard)
      hi <- if (any(binomNcode==-2)) haztemp$hindex else -1                   ## index to hx
      if (!is.null(R) && R) {
          if (!exists('simplehistoriesR')) 
              stop ("R code simplehistoriesR not available; source prwisimple.R")  
          else {
              temp <- do.call('simplehistoriesR', 
                              list (
                                  x, m, nc, cc,
                                  !is.null(MRdata$Tu), !is.null(MRdata$Tm),
                                  binomNcode, MRdata$markocc, MRdata$firstocc, pID[,x],
                                  CH, grp,
                                  array(gkhk$gk, dim=c(cc,k,m)),
                                  array(gkhk$hk, dim=c(cc,k,m)),
                                  pi.density,
                                  PIA, usge, hx, hi,
                                  maskusage))
          }
      } 
      else {
        temp <- simplehistoriescpp(
          as.integer(m),
          as.integer(nc),
          as.integer(cc),
          as.integer(grain),
          as.integer(ncores),
          as.integer(binomNcode),
          as.integer(MRdata$markocc),
          as.integer(MRdata$firstocc),
          as.double (pID[,x]),
          as.integer(CH),   
          as.integer(grp)-1L,
          as.double (gkhk$gk),     ## precomputed probability 
          as.double (gkhk$hk),     ## precomputed hazard
          as.matrix (pi.density),
          as.integer(PIA[1,,,,x]),
          as.matrix (usge),
          as.matrix (hx),                
          as.matrix (hi),      
          as.matrix (maskusage),
          as.double (telemhr),
          as.integer (telemstart))
      }
      sump   <- sump + pmixn[x,] * temp$prwi
  }
  sump 
}
#--------------------------------------------------------------------------------

integralprw1 <- function (cc0, haztemp, gkhk, pi.density, PIA0, 
  CH0, binomNcode, MRdata, grp, usge, pmixn, pID, grain, ncores) {
    nc <- dim(PIA0)[2]    ## animals
    nr <- nrow(CH0)       ## unique naive animals (1 or nc)
    m <- nrow(pi.density)
    nmix <- nrow(pmixn)
    if (length(grp)<=1) grp <- rep(1,nc)
    ngroup <- max(length(unique(grp)),1)
    sump <- numeric(nc)
    for (x in 1:nmix) {
        hx <- if (any(binomNcode==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## sum_k (hazard)
        hi <- if (any(binomNcode==-2)) haztemp$hindex else -1                   ## index to hx
        temp <- simplehistoriescpp(
          as.integer(m),
          as.integer(nr),
          as.integer(cc0),
          as.integer(grain),
          as.integer(ncores),
          as.integer(binomNcode),
          as.integer(MRdata$markocc),
          as.integer(rep(-1,nr)),                 # MRdata$firstocc  # never marked
          as.double (pID[,x]),
          as.integer(CH0),    
          as.integer(grp)-1L,              # group  
          as.double (gkhk$gk),        # precomputed probability 
          as.double (gkhk$hk),        # precomputed hazard
          as.matrix (pi.density),
          as.integer(PIA0[1,1:nr,,,x]),
          as.matrix (usge),
          as.matrix (hx),                
          as.matrix (hi),      
          as.matrix (matrix(TRUE, nrow = nr, ncol = m)),
          as.double (0),   # no telemetry
          as.integer(0)    # no telemetry
        )  
        if (nr == 1) temp$prwi <- rep(temp$prwi, nc)
        for (g in 1:ngroup) {
            ok <- as.integer(grp) == g
            sump[ok] <- sump[ok] + pmixn[x,ok] * (1-temp$prwi[ok])
        }
    }
    sump
}
#--------------------------------------------------------------------------------

expectedmu <- function (cc, haztemp, gkhk, pi.density, Nm, PIA, 
                          CH, binomNcode, MRdata, grp, usge, pmixn, pID, a0,
                          debug = FALSE) {
    nc <- nrow(CH)
    k <- nrow(usge)
    s <- ncol(usge)
    m <- nrow(pi.density)
    nmix <- nrow(pmixn)
    Tumusk <- Tmmusk <- matrix(0, k, s)
    for (x in 1:nmix) {
        hx <- if (any(binomNcode==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## lookup sum_k (hazard)
        hi <- if (any(binomNcode==-2)) haztemp$hindex else -1                   ## index to hx
        temp <- expectedmucpp(
            as.integer(nc),
            as.integer(cc),
            as.logical(!is.null(MRdata$Tu)),
            as.logical(!is.null(MRdata$Tm)),
            as.integer(MRdata$sightmodel),
            as.integer(binomNcode),
            as.integer(MRdata$markocc),
            as.double (pID[,x]),
            as.integer(grp)-1L,
            as.double (gkhk$gk),     ## precomputed probability 
            as.double (gkhk$hk),     ## precomputed hazard
            as.matrix (pi.density),
            as.matrix (Nm),
            as.integer(PIA[1,,,,x]),
            as.matrix (usge),
            as.matrix (hx),                
            as.matrix (hi),
            as.double (a0))
        Tumusk <- Tumusk + pmixn[x,1] * temp$Tumusk  ## not yet adjusted for absolute density and cell area
        Tmmusk <- Tmmusk + pmixn[x,1] * temp$Tmmusk  
    }
    list(Tumusk=Tumusk, Tmmusk=Tmmusk)
}
#--------------------------------------------------------------------------------

allhistsignal <- function (detectfn, grain, ncores, binomNcode,
                           CH, signal, grp, gk, realparval, dist2, 
                           pi.density, PIA, miscparm, maskusage,
                           pmixn, debug = FALSE) {
  nc <- nrow(CH)
  m <- nrow(pi.density)
  nmix <- nrow(pmixn)
  ngroup <- length(levels(grp))
  sump <- numeric(nc)
  for (x in 1:nmix) {
    temp <- signalhistoriescpp(
      as.integer(m),
      as.integer(nc),
      as.integer(detectfn),
      as.integer(grain),
      as.integer(ncores),
      as.integer(binomNcode),
      as.integer(CH),   
      as.matrix(signal),
      as.integer(grp)-1L,
      as.double(gk),
      as.matrix(realparval),
      as.matrix(dist2),
      as.matrix (pi.density),
      as.integer(PIA[1,1:nc,,,x]),   ## pass only PIA for x
      as.double(miscparm),
      as.matrix(maskusage))
    sump <- sump + pmixn[x,] * temp
  }
  sump
}
#--------------------------------------------------------------------------------
allhistpolygon <- function (detectfn, realparval, haztemp, hk, H, pi.density, PIA, 
                           CH, xy, binomNcode, grp, usge, mask, pmixn, maskusage,
                           grain, ncores, minprob, debug=FALSE) {
  nc <- nrow(CH)
  m <- nrow(pi.density)
  s <- ncol(usge)
  nmix <- nrow(pmixn)
  ngroup <- length(levels(grp))
  sump <- numeric(nc)
  for (x in 1:nmix) {
      hx <- if (any(binomNcode==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## lookup sum_k (hazard)
      hi <- if (any(binomNcode==-2)) haztemp$hindex else -1                   ## index to hx
      temp <- polygonhistoriescpp(
        as.integer(nc),
        as.integer(detectfn[1]),
        as.integer(grain),
        as.integer(ncores),
        as.double(minprob),          
        as.integer(binomNcode),
        as.integer(CH),   
        as.matrix(xy$xy),
        as.vector(xy$start),
        as.integer(as.numeric(grp))-1L,
        as.double(hk),
        as.double(H),
        as.matrix(realparval),
        matrix(1,nrow=s, ncol=nmix),  ## pID?
        as.matrix(mask),
        as.matrix (pi.density),
        as.integer(PIA[1,,,,x]),
        as.matrix(usge),
        as.matrix (hx),                
        as.matrix (hi),      
        as.matrix(maskusage),
        as.integer(debug)
      )
      sump <- sump + pmixn[x,] * temp
  }
  if (debug) {
    cat("SUM log(PRWI) ", sum(log(sump)), "\n")
  }
  sump
}
#--------------------------------------------------------------------------------

integralprw1poly <- function (detectfn, realparval0, haztemp, hk, H, pi.density, PIA0, 
                              CH0, binomNcode, grp, usge, mask, pmixn, maskusage,
                              grain, ncores, minprob, debug = FALSE) {
    
  nc <- dim(PIA0)[2]
  nr <- nrow(CH0)       ## unique naive animals (1 or nc)
  m <- nrow(pi.density)
  nmix <- nrow(pmixn)
  if (length(grp)<=1) grp <- rep(1,nc)
  s <- ncol(usge)
  ngroup <- length(levels(grp))
  sump <- numeric(nc)
  for (x in 1:nmix) {
      hx <- if (any(binomNcode==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## sum_k (hazard)
      hi <- if (any(binomNcode==-2)) haztemp$hindex else -1                   ## index to hx
    for (g in 1:ngroup) {
      ok <- as.integer(grp) == g
      temp <- polygonhistoriescpp(
        as.integer(nr),
        as.integer(detectfn[1]),
        as.integer(grain),
        as.integer(ncores),
        as.double(minprob),          
        as.integer(binomNcode),
        as.integer(CH0),   
        as.matrix(0L),  # empty for null history
        as.vector(0L),  # empty for null history
        as.integer(g)-1L,
        as.double(hk),
        as.double(H),
        as.matrix(realparval0),
        matrix(1,nrow=s, ncol=nmix),  ## pID?
        as.matrix(mask),
        as.matrix (pi.density),
        as.integer(PIA0[1,1:nr,,,x]),
        as.matrix(usge),
        as.matrix (hx),                
        as.matrix (hi),      
        as.matrix(maskusage),
        as.integer(debug)
      )
      if (nr == 1) temp <- rep(temp, nc)
      sump[ok] <- sump[ok] + pmixn[x,ok] * (1-temp[ok])
    }
  }
  sump
}
#--------------------------------------------------------------------------------

#######################################################################################
generalsecrloglikfn <- function (
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
  sessionLL <- function (data) {
    ## log likelihood for one session
    ## in multi-session case must get session-specific data from lists
    #---------------------------------------------------
    sessnum <- data$sessnum  # changed from argument 2021-06-22
    nc1 <- max(data$nc,1)
    PIA <- design$PIA[sessnum, 1:nc1, 1:data$s, 1:data$K, ,drop = FALSE]
    PIA0 <- design0$PIA[sessnum, 1:nc1, 1:data$s, 1:data$K, ,drop = FALSE]
    ## unmodelled beta parameters, if needed
    miscparm <- getmiscparm(details$miscparm, detectfn, beta, parindx, details$cutval)
    #---------------------------------------------------
    
    density <- getmaskpar(!CL, D, data$m, sessnum, details$unmash, 
                          attr(data$capthist, 'n.mash'))
    if (CL) {
      pi.density <- matrix(1/data$m, nrow=data$m, ncol=1)  
    }
    else {
      Dsum <- apply(density,2,sum)   ## by group
      Nm <- density * getcellsize(data$mask)
    
      if (data$MRdata$allsighting && data$MRdata$pi.mask[1] != -1) {
          pi.density <- matrix(data$MRdata$pi.mask, ncol = 1)  ## by group=column?
          if (any(Nm < nrow(data$CH)*pi.density)) {
              warning("invalid distribution for sighting at Eval ", .localstuff$iter)  # changed from stop() 2019-12-15
          }
      }
      else
          pi.density <- sweep(density, MARGIN = 2, STATS = Dsum, FUN = '/')
    }
    #---------------------------------------------------
    ## allow for scaling of detection
    Dtemp <- if (D.modelled) mean(D[,1,sessnum]) else NA
    Xrealparval <- reparameterize (realparval, detectfn, details,
                                   data$mask, data$traps, Dtemp, s)
    Xrealparval0 <- reparameterize (realparval0, detectfn, details,
                                    data$mask, data$traps, Dtemp, s)
    if (details$debug>2) browser()

    ## check valid parameter values
    if (!all(is.finite(Xrealparval))) {
      cat ('beta vector :', beta, '\n')
      warning ("extreme 'beta' in 'generalsecrloglikfn' ",
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
    
    #####################################################################
    pmixn <- getpmix (data$knownclass, PIA, Xrealparval)  ## membership prob by animal
    
    pID <- getpID(PIA, Xrealparval, data$MRdata)
    if (is.function(details$userdist)) {
      noneuc <- getmaskpar(!is.null(NE), NE, data$m, sessnum, FALSE, NULL)
      distmat2 <- getuserdist(data$traps, data$mask, details$userdist, sessnum, 
                              noneuc[,1], density[,1], miscparm, detectfn == 20)
    }
    else {
      distmat2 <- data$distmat2
    }
    ## precompute gk, hk for point detectors
    if (all(data$dettype %in% c(0,1,2,5,8,13)) || data$HPXpoly) {
        if (!is.null(details$R) && details$R) {
            if (!exists('makegkPointR')) 
                stop ("R code makegkPointR not available; source makegk.R")  
            else {
                gkhk <- do.call('makegkPointR', list(detectfn, Xrealparval, distmat2, miscparm))
            }
        }
        else {
          gkhk <- makegkPointcpp (
            as.integer(detectfn),
            as.integer(details$grain),
            as.integer(details$ncores),
            as.matrix(Xrealparval),
            as.matrix(distmat2),
            as.double(miscparm))
        }
        if (any(data$dettype == 8)) {
            gkhk <- cappedgkhkcpp (
                as.integer(nrow(Xrealparval)), 
                as.integer(nrow(data$traps)),
                as.double(attr(data$mask, "area")), 
                as.double(density[,1]), 
                as.double(gkhk$gk), 
                as.double(gkhk$hk))  
        }
    }
    ## precompute gk, hk for polygon and transect detectors
    else if (all(data$dettype %in% c(3,4,6,7))) {
      dimension <- (data$dettype[1] %in% c(3,6)) + 1   ## 1 = 1D, 2 = 2D
      # 2019-11-25 not safe to use multithreading with 2-D integration 
      # 2019-11-25 therefore using repeated 1-D integration
      convexpolygon <- is.null(details$convexpolygon) || details$convexpolygon
      gkhk <- makegkPolygoncpp (
        as.integer(detectfn), 
        as.integer(dimension), 
        as.logical(convexpolygon), 
        as.integer(details$grain), 
        as.integer(details$ncores),
        as.matrix(Xrealparval), 
        as.integer(data$cumk),
        as.matrix(data$traps), 
        as.matrix(data$mask))
      if (details$debug) {
          cat("sum(hk) ", sum(gkhk$hk), "\n")  
      }
    }
    
    ## telemetry precalculation
    if (any(data$dettype == 13)) {
        telemstart <- data$xy$start
        telemscale <- details$telemetryscale
        if (is.null(telemscale)) telemscale <- 1
        telemhr <- gethr(as.double(nrow(data$CH)),      ## or nc1?
                    as.integer(detectfn), 
                    as.double(telemstart), 
                    as.matrix(data$xy$xy), 
                    as.matrix(data$mask), 
                    as.matrix(Xrealparval), 
                    as.double(telemscale))
    }
    else {
        telemhr <- 0
        telemstart <- 0
        telemscale <- 0
    }
    #######################################################################
    ## option to estimate sighting overdispersion by simulation and exit */
    if (!is.null(details$nsim) && details$nsim > 0) {
        if (CL)
            stop("simulation for overdispersion requires full likelihood (not CL)")
        else {
          chat <- unlist(getchat (nrow(realparval0), nrow(data$CH), data$n.distrib,         ## or nc1?
            data$grp, data$usge, pmixn, pID, getcellsize(data$mask), 
            gkhk, pi.density, Dsum, PIA0, data$binomNcode, 
            data$MRdata, miscparm, details$nsim, details$grain, details$ncores))
          return (chat)         
        }
    }
    #######################################################################
    if (all(data$dettype %in% c(0,1,2,3,4,6,7,8,13))) {
        ## hazard for exclusive detectors or related
        haztemp <- gethazard (data$m, data$binomNcode, nrow(Xrealparval), gkhk$hk, PIA, data$usge)
    }
    ## model detection histories (prw) conditional on detection (pdot)
    if (data$nc == 0) {
        prw <- 1  ## simple if no animals detected
    }
    else {
        if (all(data$dettype %in% c(0,1,2,8,13)) || data$HPXpoly) {
            prw <- allhistsimple (nrow(Xrealparval), haztemp, gkhk, pi.density, PIA, 
                                  data$CH, data$binomNcode, data$MRdata, data$grp, data$usge, pmixn, 
                                  pID, data$maskusage, 
                                  telemhr, telemstart, telemscale,
                                  details$grain, details$ncores, details$R)
        }
        else if (all(data$dettype == 5)) {
            prw <- allhistsignal (detectfn, details$grain, details$ncores, data$binomNcode, data$CH, data$signal$signal,
                                  data$grp, gkhk$gk, Xrealparval, distmat2, pi.density, PIA, 
                                  miscparm, data$maskusage, pmixn)
        }
        else if (all(data$dettype %in% c(3,4,6,7))) {
            prw <- allhistpolygon (detectfn, Xrealparval, haztemp, gkhk$hk, gkhk$H, pi.density, PIA, 
                                   data$CH, data$xy, data$binomNcode, data$grp, data$usge, data$mask,
                                   pmixn, data$maskusage, details$grain, details$ncores, details$minprob,
                                   debug = details$debug>3)
        }
        else {
            stop ("this detector type, or mixed detector types, not available yet in secr 4.6")
        }
    }    
        ## polygon types
    if (all(data$dettype %in% c(3,4,6,7)) && !data$HPXpoly) {
        if (learnedresponse) {   ## overwrite gk,hk with model for naive animal
            gkhk <- makegkPolygoncpp (
              as.integer(detectfn), 
              as.integer(details$grain),
              as.integer(details$ncores),
              as.matrix(Xrealparval0),
              as.integer(data$cumk),
              as.matrix(data$traps),
              as.matrix(data$mask))
            if (all(data$dettype %in% c(3,4))) {
                ## hazard for exclusive detectors or related bug fix 2020-04-24
                haztemp <- gethazard (data$m, data$binomNcode, nrow(Xrealparval0), gkhk$hk, PIA0, data$usge)
            }
        }
        pdot <- integralprw1poly (detectfn, Xrealparval0, haztemp, gkhk$hk, gkhk$H, pi.density, PIA0, 
                                  data$CH0, data$binomNcode, data$grp, data$usge, data$mask,
                                  pmixn, data$maskusage, details$grain, details$ncores, details$minprob, 
          debug = details$debug>3)
    }
    ## point types
    else {
        if (learnedresponse) {   ## overwrite gk,hk with model for naive animal
            if (!is.null(details$R) && details$R) {   # inserted 2020-04-24 to use R code here, too
                if (!exists('makegkPointR')) 
                    stop ("R code makegkPointR not available; source makegk.R")  
                else {
                    gkhk <- do.call('makegkPointR', list(detectfn, Xrealparval0, distmat2, miscparm))
                }
            }
            else gkhk <- makegkPointcpp (
                as.integer(detectfn), 
                as.integer(details$grain),
                as.integer(details$ncores),
                as.matrix(Xrealparval0), 
                as.matrix(distmat2), 
                as.double(miscparm))
            ## no capped adjustment as learned response not compatible
            
            if (all(data$dettype %in% c(0,8))) {
                ## hazard for exclusive detectors or related bug fix 2020-04-24
                haztemp <- gethazard (data$m, data$binomNcode, nrow(Xrealparval0), gkhk$hk, PIA0, data$usge)
            }
            
        }
        pdot <- integralprw1 (nrow(Xrealparval0), haztemp, gkhk, pi.density, PIA0, 
                              data$CH0, data$binomNcode, data$MRdata, data$grp, data$usge, pmixn, 
                            pID, details$grain, details$ncores)
    }
    
    ngroup <- max(length(levels(data$grp)),1)
    comp <- matrix(0, nrow = 6, ncol = ngroup)
    for (g in 1:ngroup) {
      ok <- as.integer(data$grp) == g
      #----------------------------------------------------------------------
      
      comp[1,g] <- if (any(is.na(prw)) || any(prw<=0)) NA else sum(log(prw[ok]))
      
      #----------------------------------------------------------------------
      ## Adjust for undetected animals unless data includes all-zero histories
      ## (the case for allsighting data when knownmarks = TRUE).
      if (!data$MRdata$sightmodel==5 && !all(data$dettype==13)) {
          comp[2,g] <- if (any(is.na(pdot)) || any(pdot<=0)) NA else -sum(log(pdot[ok]))
      }
      
      #----------------------------------------------------------------------
      
      if (!CL && !data$MRdata$allsighting) {
          ng <- sum(ok)
          if (any(data$dettype==13))
              nonzero <- sum(apply(data$CH[,data$dettype!=13,,drop=FALSE] != 0,1,sum)[ok]>0)
          else
              nonzero <- ng
          N <- sum(Nm[,g])
          if (ng == 0) {
              meanpdot <- pdot
          }
          else {
              meanpdot <- ng / sum(1/pdot[ok])
          }
          ## 2023-09-22
          if (data$n.distrib == 1 && .localstuff$iter == 0 && nonzero>N) {
              warning("distribution = 'binomial' ",
                      "but number detected n (", nonzero, 
                      ") exceeds initial value of N (", round(N,1), ")")
          }
              
          comp[3,g] <- if (is.na(meanpdot) || (meanpdot <= 0)) NA 
              else switch (data$n.distrib+1,
                               dpois(nonzero, N * meanpdot, log = TRUE),
                               lnbinomial (nonzero, N, meanpdot),
                               NA)
      }
      #----------------------------------------------------------------------
      # adjustment for mixture probabilities when class known
      known <- sum(data$knownclass[ok]>1)
      if (details$nmix>1 && known>0) {
          nb <- details$nmix + 1
          nm <- tabulate(data$knownclass[ok], nbins = nb)
          pmix <- attr(pmixn, 'pmix')
          
          # for (x in 1:details$nmix) {
          #     # need group-specific pmix
          #     # comp[4,g] <- comp[4,g] + nm[x+1] * log(pmix[x]) 
          # }
          
          ## 2022-01-16 bug fix
          # firstx <- match ((1:details$nmix)+1, data$knownclass[ok])
          # tempsum <- sum(pdot[ok][firstx] * pmix)
          # comp[4,g] <- sum(nm[-1] * log(pdot[ok][firstx] * pmix / tempsum))
          
          ## 2022-10-25 bug fix
          firstx <- match ((1:details$nmix)+1, data$knownclass)
          pdpmix <- pdot[firstx] * pmix
          pdpmix <- pdpmix[!is.na(pdpmix)]
          comp[4,1] <- sum(nm[-1] * log(pdpmix / sum(pdpmix)))
          
      }
   
      #----------------------------------------------------------------------
      # sightings
      sightingocc <- data$MRdata$markocc < 1
      if (any(sightingocc)) {
          Nm <- density * getcellsize(data$mask)
          tmp <- expectedmu (nrow(Xrealparval), haztemp, gkhk, pi.density, Nm, PIA, 
                             data$CH, data$binomNcode, data$MRdata, data$grp, data$usge, pmixn, 
                             pID, pdot[1])
          Tumusk <- tmp$Tumusk ## * sum(density[,g]) * getcellsize(data$mask)
          Tmmusk <- tmp$Tmmusk ## * sum(density[,g]) * getcellsize(data$mask)
          if (!is.null(data$MRdata$Tu) && !is.null(Tumusk)) {
              Tu <- data$MRdata$Tu
              Tulik <- Tsightinglikcpp (Tu, data$MRdata$markocc, data$binomNcode,
                                        data$usge, Tumusk, details$debug)
              if (Tulik$resultcode != 0) 
                  comp[5,1] <- NA
              else
                  comp[5,1] <- Tulik$Tlik/details$chat[1] 
          }
          if (!is.null(data$MRdata$Tm) && !is.null(Tmmusk)) {
              Tm <- data$MRdata$Tm
              Tmlik <- Tsightinglikcpp (Tm, data$MRdata$markocc, data$binomNcode,
                                        data$usge, Tmmusk, details$debug)
              if (Tmlik$resultcode != 0)
                  comp[6,1] <- NA
              else
                  comp[6,1] <- Tmlik$Tlik/details$chat[2]
          }
      }
      #----------------------------------------------------------------------
    }   ## end loop over groups
    if (details$debug>=1) {
        ## display likelihood components summed over groups, and logmultinomial constant
        comp <- apply(comp,1,sum)
        cat(comp[1], comp[2], comp[3], comp[4], comp[5], comp[6], data$logmult, '\n')
    }
    sum(comp) + data$logmult
  
  } ## end sessionLL
  
  ######################################################################################
  ## Main line of generalsecrloglikfn
  ######################################################################################
  if (details$debug>4) browser()
  nsession <- length(sessionlevels)
  #--------------------------------------------------------------------
  # Fixed beta
  beta <- fullbeta(beta, details$fixedbeta)
  #--------------------------------------------------------------------
  # Detection parameters
  detparindx <- parindx[!(names(parindx) %in% c('D', 'noneuc'))]
  detlink <- link[!(names(link) %in% c('D', 'noneuc'))]
  realparval  <- makerealparameters (design, beta, detparindx, detlink, fixed)
  realparval0 <- makerealparameters (design0, beta, detparindx, detlink, fixed)
  #--------------------------------------------------------------------
  sessmask <- lapply(data, '[[', 'mask')
  grplevels <- unique(unlist(lapply(data, function(x) levels(x$grp))))
  #---------------------------------
  # Density
  D.modelled <- !CL & is.null(fixed$D)
  if (!CL ) {
      D <- getD (designD, beta, sessmask, parindx, link, fixed,
               grplevels, sessionlevels, parameter = 'D')
  }
  #--------------------------------------------------------------------
  # Non-Euclidean distance parameter
  NE <- getD (designNE, beta, sessmask, parindx, link, fixed,
              grplevels, sessionlevels, parameter = 'noneuc')

  #--------------------------------------------------------------------
  # Two types of call
  # (i) overdispersion of sightings simulations only
  if (details$nsim > 0) {   
    ## chat <- mapply (sessionLL, data, 1:nsession, SIMPLIFY = FALSE)
    chat <- mapply (sessionLL, data, SIMPLIFY = FALSE)
    chatmat <- matrix(unlist(chat), ncol = 3, byrow = TRUE)
    dimnames(chatmat) <- list(session = 1:nsession, chat = c('Tu', 'Tm','Tn'))
    return(chatmat)
  }
  #--------------------------------------------------------------------
  # (ii) typical likelihood evaluation
  else {
    loglik <- sum(sapply (data, sessionLL)) 
    .localstuff$iter <- .localstuff$iter + 1  
      if (details$trace) {
          fixedbeta <- details$fixedbeta
          if (!is.null(fixedbeta))
              beta <- beta[is.na(fixedbeta)]
          cat(format(.localstuff$iter, width=4),
              formatC(round(loglik,dig), format='f', digits=dig, width=10),
              formatC(beta, format='f', digits=dig+1, width=betaw),
              '\n')
          flush.console()
      }
      loglik <- ifelse(is.finite(loglik), loglik, -1e10)
      ifelse (neglik, -loglik, loglik)
  }
}  ## end of generalsecrloglikfn
############################################################################################

