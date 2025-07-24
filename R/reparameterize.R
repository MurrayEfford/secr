###############################################################################
## package 'secr'
## reparameterize.R

reparameterize.sigmak <- function (realparval, D, linear) {
  ## D,sigmak parameterisation 2014-03-12
  ## 2021-06-20 included 'd' parameter 
  ## vector D must match rows in realparval
  realnames <- dimnames(realparval)[[2]]
  sigmakindex <- match('sigmak', realnames)
  cindex <- match('c', realnames)
  dindex <- match('d', realnames)
  if (is.na(cindex))
    cval <- 0
  else
    cval <- realparval[,cindex]
  if (is.na(dindex))
    dval <- 0
  else
    dval <- realparval[,dindex]
  if (!('sigmak' %in% realnames))
    stop ("'param = 4:6 ' requires 'sigmak' in model")
  if (linear)
    realparval[,sigmakindex] <- realparval[,sigmakindex] / (D+dval) + cval
  else
    realparval[,sigmakindex] <- realparval[,sigmakindex] / (D+dval)^0.5 + cval
  dimnames(realparval)[[2]][sigmakindex] <- 'sigma'
  realparval <- realparval[, -c(cindex, dindex), drop = FALSE]   ## added 2014-08-18, modified 2021-06-20
  realparval
}
###############################################################################

reparameterize.esa <- function (realparval, mask, traps, detectfn, nocc) {
  
  ## esa, sigma parameterisation 2013-07-01
  ## could whole fn be coded in C?
  g0fromesa <- function (a, sigma, z, lower = 0, upper = 1) {
    fx <- function(g0) {
      ## pdot accounts for 'usage'
      ## pdot selects appropriate g0/lambda0 according to detectfn
      ## 2021-05-21 remove ncores = 1
      (sum(pdot(mask, traps, detectfn = detectfn,
        detectpar = list(g0 = g0, lambda0 = g0, sigma = sigma, z = z),
        noccasions = nocc)) * cell) - a
    }
    tmp <- try(uniroot(fx, lower=lower, upper=upper), silent = TRUE)
    ## debug if (inherits(tmp, 'try-error')) print(c(fx(lower),fx(upper)))
    if (inherits(tmp, 'try-error')) 0.0001 else tmp$root
  }
  cell <- secr_getcellsize(mask)
  if (inherits(mask, 'linearmask'))
    stop ('esa parameterization not available for linear masks')
  dettype <- secr_detectorcode(traps)
  if (!all(dettype %in%  c(0,1,2,5,8,13))) {
      stop ("esa parameterization is available only for point detectors")
  }
  realnames <- dimnames(realparval)[[2]]
  sigmaindex <- match('sigma', realnames)
  esaindex <- match('esa', realnames)
  ndetectpar <- length(secr_parnames(detectfn))
  z <- ifelse (ndetectpar == 3, realparval[, match('z', realnames)], 1)
  if (is.na(esaindex) | is.na(sigmaindex))
    stop ("'param = 2' requires both 'esa' and 'sigma' in model")
  realparval[,esaindex]  <- unlist(mapply(g0fromesa,
                                          realparval[,esaindex],  ## a
                                          realparval[,sigmaindex],
                                          z
  ))
  realparval
}
###############################################################################

reparameterize.a0 <- function (realparval, detectfn, linear) {
  ## a0, sigma parameterisation 2013-07-19
  realnames <- dimnames(realparval)[[2]]
  sigmaindex <- match('sigma', realnames)
  a0index <- match('a0', realnames)
  if (! all (c('a0','sigma') %in% realnames))
    stop ("'param = 3 or 5 ' requires both 'a0' and 'sigma' in model")
  if (!(detectfn %in% c(0:8, 14:19)))
    stop ('invalid combination of param = 3 or 5 and detectfn')
  
  if (linear)
    lambda0 <- realparval[,a0index] / realparval[,sigmaindex] * 1000
  else
    lambda0 <- realparval[,a0index] / 2 / pi / realparval[,sigmaindex]^2 * 10000
  realparval[,a0index] <- if (detectfn %in% 0:8) 1-exp(-lambda0) else lambda0
  dimnames(realparval)[[2]][a0index] <- if (detectfn<9) 'g0' else 'lambda0'
  realparval
}
###############################################################################

secr_reparameterize <- function (realparval, detectfn, details, mask, traps, D, s) {
  
  ##----------------------------------------------
  ## allow for scaling of detection in one session
  
  ## D is scalar density or NA
  ## s is number of occasions
  
  linear <- inherits(mask, 'linearmask')
  if (details$param %in% 4:6) {
    ## does not allow varying density surface
    ## cf scaled.detection()
    realparval <- reparameterize.sigmak (realparval, D, linear)
  }
  
  if (details$param %in% c(2,6)) {
    realparval <- reparameterize.esa (realparval, mask, traps, detectfn, s)
  }
  else if (details$param %in% c(3,5)) {
    realparval <- reparameterize.a0 (realparval, detectfn, linear)
  }
  
  if (all(detector(traps) == 'telemetry'))
    realparval[,'lambda0'] <- 1.0
  
  realparval
}
###############################################################################

