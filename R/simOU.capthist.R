# Do not need mvtnorm when x,y uncorrelated

simOU <- function(xy, tau, sigma, noccasions, start = NULL){
    # xy: long-term mean (mu) for both dimensions
    # tau: time constant (1/beta), determines mean reversion speed
    # sigma: sd of stationary distribution (not SDE diffusion parameter)
    # noccasions: number of simulation steps
    # start: initial position (optional, defaults to stationary distribution)
    
    if (!is.null(dim(xy)) && nrow(xy)>1) {
        apply(xy, 1, simOU, simplify = FALSE,
              tau = tau, sigma = sigma, noccasions = noccasions, start = start)
    }
    else {
        xy <- as.numeric(xy)   # dodge issue with dataframe
        
        # If start is not provided, initialize from the stationary distribution
        if (is.null(start)) {
            start <- rnorm(2, mean = xy, sd = sigma)
        }
        
        # Calculate the standard deviation for the discrete time step
        # assuming Delta t = 1 for each step in this discrete simulation
        beta <- 1/tau
        step_sd <- sigma * sqrt(1 - exp(-2*beta))
        
        out <- matrix(0, nrow = noccasions, ncol = 2)
        out[1, ] <- start
        for (i in 2:noccasions){
            # Mean calculation for the current step (from analytical solution)
            current_mean <- (1 - exp(-beta)) * xy + exp(-beta) * out[i - 1, ]
            # Draw 2 independent normal random variates for the x and y dimensions
            out[i, ] <- rnorm(2, mean = current_mean, sd = step_sd)
        }
        out
    }
}

# verify 2025-06-13

# sigma is sd of stationary distribution, not SDE diffusion coefficient
# hence absorbs 1/(2 beta) in variance

# nrepl <- 1000
# taulevels <- c(0.001, 0.1,1,10,100)
# nocclevels <-  c(10,20,100,1000,5000)
# out <- array(dim=c(5,5,nrepl), dimnames=list(tau=taulevels, nocc=nocclevels, NULL))
# for (k in 1:nrepl) {
#     for (i in 1:5) {
#         for (j in 1:5) {
#             locs <- simOU(c(0,0), tau = taulevels[i], sigma = 1,
#                           noccasions = nocclevels[j])
#             out[i,j,k] <- mean(apply(locs,2,var))^0.5  # RMS
#         }
#     }
# }
# 
# round(apply(out,1:2,mean),4)

#         nocc
# tau         10     20    100   1000   5000
#   0.001 0.9870 0.9925 0.9975 1.0005 0.9995
#   0.1   0.9894 0.9908 0.9976 0.9998 1.0003
#   1     0.9323 0.9641 0.9922 0.9998 0.9999
#   10    0.5159 0.6545 0.9019 0.9893 0.9977
#   100   0.1788 0.2481 0.5044 0.8963 0.9797

simOU.capthist <- function (
        traps,
        popn,
        detectpar,     # list of epsilon, sigma, tau
        noccasions,    # effective "duration"
        seed     = NULL,
        savepopn = FALSE,
        savepath = FALSE,
        ...)
{
    ##################
    ## set random seed
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
    
    captfn <- function (xy) edist(xy,traps) <= detectpar$epsilon   
    N <- nrow(popn)
    locs <- apply(popn, 1, simOU, detectpar$tau, detectpar$sigma, noccasions, simplify = FALSE)
    capt <- lapply(locs, captfn)
    capt <- do.call(rbind, capt)  
    ch <- array(capt, dim = c(noccasions, N, ncol(capt)), 
                dimnames = list(1:noccasions,rownames(popn),rownames(traps)))
    ch <- aperm(ch, c(2,1,3))
    ch <- ch[apply(ch, 1, sum) > 0,,, drop = FALSE]   # drop null histories
    class(ch) <- 'capthist'
    traps(ch) <- traps
    # cast as required detector type
    if (detector(traps)[1] != "multi" || length(list(...))>0) {
        ch <- reduce(ch, outputdetector = detector(traps)[1], dropunused = FALSE, ...)
    }
    attr(ch, 'seed')      <- RNGstate      ## save random seed
    attr(ch, 'detectpar') <- detectpar
    if (savepopn) attr(ch, 'popn') <- popn
    if (savepath) attr(ch, 'path') <- locs
    if (verify) verify(ch)
    ch
}
