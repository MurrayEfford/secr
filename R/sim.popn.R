###############################################################################
## package 'secr'
## sim.popn.R
## simulate spatially distributed population
## transferred from methods.R 2010-06-08
## last changed
## 2010-06-14  added Ndist = 'specified', using D to carry N
## 2011-03-27  conditional attachment of rownames (correct bug when N = 0)
## 2011-04-06  coastal beta option
## 2011-10-20  poly argument
## 2012-04-10  IHP D may be covariate; more checks
## 2012-04-10  MRC
## 2013-11-23 IHP for Ndist fixed
## 2014-04-18 session-specific density
## 2014-09-03 model2D = "linear" option, and some recoding of IHP
## 2014-12-29 buffertypes concave, convex
## 2015-01-13 debugged non-rectangular buffers; see sim.popn.test.R in testing
## 2015-02-18 multisession sim.popn updated for Nbuffer
## 2016-09-21 model2D = "even" option
## 2017-06-07 missing D becomes NULL
## 2018-02-21 multinomial recruitment
## 2018-06-26 stop.at.edge
## 2020-11-09 radialexp, radialnorm undocumented movement options
## 2020-12-11 "normalize" option for edgemethod; internal 'disperse' function
## 2021-02-12 "clipandreplace" edgemethod option; clip can be used with constantN
## 2021-02-13 corrected bug in "normalize" option for edgemethod
## 2021-06-01 BVE, BVN, BVT aliases for kernel names
## 2021-07-12 movement simulations all polar coordinates r,theta:
## 2021-07-12 BVE simulation uses rgamma
## 2021-07-12 BVN simulation uses rweibull
## 2021-07-12 BVT simulation uses inverse distribution function
## 2021-09-08 "truncate" as synonym of "normalize"
## 2022-06-06 IHP safe for multicolumn D df)
## 2023-05-30 IHP rmultinom handles boundary N = 0
## 2023-08-19 model2D = "rLGCP"
## 2023-08-21 model2D = "rThomas"
## 2023-10-10 details$clone == 'constant' only fixes n offspr, scale still applies
## 2023-11-04 save parents rThomas
###############################################################################

toroidal.wrap <- function (pop) {
    bb <- attr(pop, 'boundingbox')
    xmin <- min(bb$x)
    xmax <- max(bb$x)
    ymin <- min(bb$y)
    ymax <- max(bb$y)
    xrange <- xmax-xmin
    yrange <- ymax-ymin
    remainder <- function (x,y) x - (x %/% y * y)
    pop$x <- ifelse (pop$x>xmax, xmin + remainder (pop$x-xmax, xrange), pop$x)
    pop$x <- ifelse (pop$x<xmin, xmax - remainder (xmin - pop$x, xrange), pop$x)
    pop$y <- ifelse (pop$y>ymax, ymin + remainder (pop$y-ymax, yrange), pop$y)
    pop$y <- ifelse (pop$y<ymin, ymax - remainder (ymin - pop$y, yrange), pop$y)
    pop
}
# drop.outside <- function (pop) {
#     subset(pop, inside(pop))
# }

## Based on 
## https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
## Returns c(NA,NA) if the lines do not intersect, otherwise x,y vector of intersection point
get.line.intersection <- function (p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, p3.x, p3.y)
{
    s1.x <- p1.x - p0.x
    s1.y <- p1.y - p0.y
    s2.x <- p3.x - p2.x
    s2.y <- p3.y - p2.y
    
    s <- (-s1.y * (p0.x - p2.x) + s1.x * (p0.y - p2.y)) / (-s2.x * s1.y + s1.x * s2.y)
    t <- ( s2.x * (p0.y - p2.y) - s2.y * (p0.x - p2.x)) / (-s2.x * s1.y + s1.x * s2.y)
    
    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        ## Collision detected
        i.x <- p0.x + (t * s1.x)
        i.y <- p0.y + (t * s1.y)
        return (c(i.x, i.y))
    }
    return (c(NA,NA))  ## No collision
}

## vectorized
get.line.intersection.v <- function (p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, p3.x, p3.y)
{
    s1.x <- p1.x - p0.x
    s1.y <- p1.y - p0.y
    s2.x <- p3.x - p2.x
    s2.y <- p3.y - p2.y
    s <- (-s1.y * (p0.x - p2.x) + s1.x * (p0.y - p2.y)) / (-s2.x * s1.y + s1.x * s2.y)
    t <- ( s2.x * (p0.y - p2.y) - s2.y * (p0.x - p2.x)) / (-s2.x * s1.y + s1.x * s2.y)
    out <- cbind(i.x = p0.x + (t * s1.x), i.y = p0.y + (t * s1.y))
    OK <- (s >= 0 & s <= 1 & t >= 0 & t <= 1)
    out[!OK,] <- NA
    out
}

stop.at.edge <- function (old, pop, tol = 1e-6) {
    ## animals that stop at the edge are shifted tol m towards centre 
    ## so they don't get 'stuck'
    update<- function(pop1, pop2) {
        stopped <- !is.na(pop1[,1])
        if (any(stopped))
            pop2[stopped, ] <- pop1[stopped,]
        pop2
    }
    bb <- attr(pop, 'boundingbox')
    xmin <- min(bb$x)
    xmax <- max(bb$x)
    ymin <- min(bb$y)
    ymax <- max(bb$y)
    ## warning("get.line.intersection.v not working 2018-06-26")
    left <- get.line.intersection.v(pop$x, pop$y, old$x, old$y, xmin, ymin, xmin, ymax)
    right <- get.line.intersection.v(pop$x, pop$y, old$x, old$y, xmax, ymin, xmax, ymax)
    top <- get.line.intersection.v(pop$x, pop$y, old$x, old$y, xmin, ymax, xmax, ymax)
    bottom <- get.line.intersection.v(pop$x, pop$y, old$x, old$y, xmin, ymin, xmax, ymin)
    left[,1] <- left[,1] + tol
    right[,1] <- right[,1] - tol
    top[,2] <- top[,2] - tol
    bottom[,2] <- bottom[,2] + tol
    pop <- update(left, pop)
    pop <- update(right, pop)
    pop <- update(top, pop)
    pop <- update(bottom, pop)
    pop
}

reflect <- function (pop) {
    bb <- attr(pop, 'boundingbox')
    xmin <- min(bb$x)
    xmax <- max(bb$x)
    ymin <- min(bb$y)
    ymax <- max(bb$y)
    OK <- (pop$x>=xmin) & (pop$x<=xmax) & (pop$y>=ymin) & (pop$y<=ymax)
    if (!all(OK)) {
        while (any(pop$x < xmin)) pop$x <- ifelse(pop$x<xmin, 2*xmin-pop$x, pop$x)
        while (any(pop$x > xmax)) pop$x <- ifelse(pop$x>xmax, 2*xmax-pop$x, pop$x)
        while (any(pop$y < ymin)) pop$y <- ifelse(pop$y<ymin, 2*ymin-pop$y, pop$y)
        while (any(pop$y > ymax)) pop$y <- ifelse(pop$y>ymax, 2*ymax-pop$y, pop$y)
    }
    pop
}

tile <- function (popn, method = "reflect") {
    bbox <- attr(popn, 'boundingbox')
    if (method== "reflect") {
        p2 <- rbind(popn, flip(popn,lr=min(bbox$x)), flip(popn,lr=max(bbox$x)))
        rbind(p2, flip(p2,tb=min(bbox$y)), flip(p2,tb=max(bbox$y)))
    }
    else if (method == "copy") {
        ht <- max(bbox$y) - min(bbox$y)
        wd <- max(bbox$x) - min(bbox$x)
        p2 <- rbind(popn, shift(popn,c(-wd,0)), shift(popn, c(wd,0)))
        rbind(p2, shift(p2,c(0,-ht)), shift(p2, c(0,ht)))
    }
    else
        stop ("unrecognised method")
}

# post-dispersal locations
# see main function sim.popn() for application of edge methods
disperse <- function (newpopn, turnoverpar, t, core, disp) {
    nsurv <- sum(disp)
    move.a <- turnoverpar$move.a[t]
    move.b <- turnoverpar$move.b[t]
    
    if (is.function(turnoverpar$movemodel)) {
        f <- formals(turnoverpar$movemodel)
        if (length(f)==2) {
            newpopn[disp,] <- newpopn[disp,] + turnoverpar$movemodel(nsurv, move.a)
        }
        else if (length(f)==3) {
            newpopn[disp,] <- newpopn[disp,] + turnoverpar$movemodel(nsurv, move.a, move.b)
        }
        else stop ("invalid movement function")
    }
    else {
        if (turnoverpar$movemodel %in% c('normal', 'BVN')) {
            r <- rweibull(nsurv, shape = 2, scale = sqrt(2) * move.a)
        }
        else if (turnoverpar$movemodel %in% c('exponential', 'BVE')) {
            r <- rgamma(nsurv, shape = 2, scale = move.a)   
        }
        else if (turnoverpar$movemodel %in% c('t2D', 'BVT')) {
            Finv <- function (u, alpha, beta) sqrt((u^(-1/beta) - 1) * alpha^2)
            r <- Finv(runif(nsurv), move.a, move.b)
        }
        else if (turnoverpar$movemodel %in% c('frE','RDE')) {
            r <- rexp(nsurv, rate = 1/move.a)  
        }
        else if (turnoverpar$movemodel %in% c('frG','RDG')) {
            r <- rgamma(nsurv, shape = move.b, rate = 1/move.a)   
        }
        else if (turnoverpar$movemodel %in% c('frL','RDL')) {
            mu <- log(move.a)
            s <- sqrt(log(1 + 1/move.b))
            r <- rlnorm(nsurv, meanlog = mu, sdlog = s)  
        }
        else if (turnoverpar$movemodel %in% c('uniform','UNI')) {
            r <- runif(nsurv)^0.5 * move.a
        }
        else if (turnoverpar$movemodel %in% c('frZ')) {
            r <- ifelse(runif(nsurv)<move.b, 0, turnoverpar$d0 + rexp(nsurv, rate = 1/move.a))
        }
        else if (turnoverpar$movemodel == 'radialexp') {
            centre <- apply(core,2,mean)
            dxy <- sweep(newpopn[disp,], STATS = centre, FUN = "-", MARGIN = 2)
            theta <- atan2(dxy[,2],dxy[,1]) 
            r <- rgamma(nsurv, shape = 2, scale = move.a)
        }
        else if (turnoverpar$movemodel == 'radialnorm') {
            centre <- apply(core,2,mean)
            dxy <- sweep(newpopn[disp,], STATS = centre, FUN = "-", MARGIN = 2)
            theta <- atan2(dxy[,2],dxy[,1])
            r <- abs(rnorm (nsurv, mean = 0, sd = move.a))
        }
        else {
            r <- 0
            warning("unsupported movement model ", turnoverpar$movemodel, " ignored")
        }
        
        theta <- runif(nsurv, 0, 2 * pi)
      
        # 2021-11-23,28
        candidates <- newpopn[disp,] + r * cbind(cos(theta), sin(theta))
        if (inherits(core,'mask') && !is.null(covariates(core)$settle)) {
            tmp <- addCovariates(candidates, core, 'settle') 
            # assign 'outside mask' location to animals yet-to-settle
            candidates[runif(nsurv) > covariates(tmp)$settle,] <- -Inf
        }
        
        newpopn[disp,] <- candidates
    }
    newpopn
}

sim.popn <- function (D, core, buffer = 100, model2D = c("poisson", 
    "cluster", "IHP", "coastal", "hills", "linear", "even", "rLGCP", "rThomas"), 
    buffertype = c("rect", "concave", "convex"), poly = NULL,
    covariates = list(sex = c(M = 0.5,F = 0.5)), number.from = 1, Ndist
    = c('poisson','fixed','specified'), nsessions = 1, details = NULL,
    seed = NULL, keep.mask = model2D %in% c('IHP','linear'), Nbuffer = NULL,
    age = FALSE,
    ...)  {
    inside <- function(pop) {
        if (model2D != 'IHP') {
            bb <- attr(pop, 'boundingbox')
            xmin <- min(bb$x)
            xmax <- max(bb$x)
            ymin <- min(bb$y)
            ymax <- max(bb$y)
            (pop$x>=xmin) & (pop$x<=xmax) & (pop$y>=ymin) & (pop$y<=ymax)
        }
        else {
            # assuming core is a mask
            pointsInPolygon(pop, core)
        }
    }
    
    if (missing(D)) D <- NULL
    model2D <- match.arg(model2D)
    Ndist <- match.arg(Ndist)
    buffertype <- match.arg(buffertype)
    if (buffertype %in% c('convex','concave') & (model2D != 'poisson'))
        stop ("buffertype incompatible with model2D")
    if (model2D == 'even' && Ndist != 'fixed') {
        warning ('Ndist is coerced to "fixed" when model2D even')
        Ndist <- 'fixed'
    }
    if (model2D %in% c("rLGCP", "rThomas") && Ndist == 'fixed') {
        warning ('Ndist is coerced to "poisson" when model2D rLGCP, rThomas')
        Ndist <- 'poisson'
    }
    lastnumber <- number.from-1
    if (nsessions > 1) {
        discrete <- function(x) {
            fr <- x-trunc(x)
            sample (c(trunc(x), trunc(x)+1), size=1, prob=c(1-fr, fr))
        }
        session.popn <- function (s, D=NULL, Nbuffer=NULL, Ndist) {
            ## independent population
            if (s > 1) seed <- NULL   ## 2015-02-18
            if (!is.null(Nbuffer) && is.na(Nbuffer)) {
                Nbuffer <- NULL
            }
            ## sim.popn (D[1], core, buffer, model2D, buffertype, poly,
            if (ms(core)) core <- core[[s]]
            sim.popn (D, core, buffer, model2D, buffertype, poly,
                covariates, number.from = lastnumber+1, Ndist, nsessions = 1, details, seed,
                keep.mask, Nbuffer[1])
        }
        turnover <- function (oldpopn, t) {
            ## project existing population
            ## assume normal movement kernel
            ## assume lambda lacks process variance
            ## ideally lambda lognormal
            ## need 'wrap' option for toroidal wrapping of 'rect' locations
            # newstart <- max(as.numeric(rownames(oldpopn))) + 1
            if (turnoverpar$survmodel=='binomial') {
                survive <- sample (c(FALSE, TRUE), nrow(oldpopn), replace = TRUE,
                                   c(1-turnoverpar$phi[t],turnoverpar$phi[t]))
                nsurv <- sum(survive)
            }
            else {   ## assume 'discrete'
                nsurv <- discrete (turnoverpar$phi[t] * nrow(oldpopn))
                survive <- sample (nrow(oldpopn), replace = FALSE, size = nsurv)
                survive <- sort(survive)   ## numeric indices
            }
            
            if (!is.function(turnoverpar$movemodel) &&
                    turnoverpar$movemodel %in% c('uncorrelated','IND')) {
                newpopn <- sim.popn(D = D, core = core, buffer = buffer,
                    model2D = model2D, buffertype = buffertype, poly = poly,
                    covariates = covariates, Ndist = 'specified', Nbuffer = nsurv,
                    nsessions = 1, details = details)
                row.names(newpopn) <- row.names(oldpopn)[survive]
                
                ## 2021-10-15 INDzi
                if (turnoverpar$zeroinflated) {
                    stayhome <- runif(nrow(newpopn)) < turnoverpar$move.a[t]
                    newpopn[stayhome,] <- oldpopn[survive,][stayhome,]
                }
                    
            }
            else {
                newpopn <- subset(oldpopn, subset = survive)
                if (is.function(turnoverpar$movemodel) || 
                        !turnoverpar$movemodel %in% c('static','uncorrelated', 'IND')) {
                    if (turnoverpar$move.a[t] > 0) {      ## condition revived 2020-11-09
                        oldposition <- newpopn  ## remember starting position
                        dispersing <- rep(TRUE, nsurv)
                        #-------------------------------------------------------
                        newpopn <- disperse(newpopn, turnoverpar, t, core, dispersing)
                        #-------------------------------------------------------
                        if (turnoverpar$zeroinflated) {
                            stayhome <- runif(nrow(newpopn)) < turnoverpar$move.b[t]
                            newpopn[stayhome,] <- oldposition[stayhome,]
                        }
                        #-------------------------------------------------------
                        
                        if (turnoverpar$edgemethod == "wrap") {
                            newpopn <- toroidal.wrap(newpopn)
                        }
                        else if (turnoverpar$edgemethod %in% c("clip", "clipandreplace")) {
                            emigrant <- !inside(newpopn)
                            if (turnoverpar$edgemethod == "clip") {
                                newpopn <- subset(newpopn, !emigrant)
                                # new 2021-02-12 to allow constantN with "clip"
                                nsurv <- nrow(newpopn)
                            }
                            # new 2021-02-12 replace each emigrant with an immigrant at the same spot
                            if (turnoverpar$edgemethod == "clipandreplace") {
                                nonempopn <- subset(newpopn, !emigrant)
                                immpopn <- subset(oldposition, emigrant)
                                nimm <- nrow(immpopn)
                                rownames(immpopn) <- lastnumber + (1:nimm)
                                newpopn <- rbind(nonempopn, immpopn, renumber = FALSE)
                                lastnumber <<- lastnumber + nimm
                            }
                        }
                        else if (turnoverpar$edgemethod == "stop") {
                            newpopn <- stop.at.edge(oldposition, newpopn)
                        }
                        else if (turnoverpar$edgemethod == "reflect") {
                            newpopn <- reflect(newpopn)
                        }
                        else if (turnoverpar$edgemethod %in% c("normalize", "truncate")) {
                            dispersing <- !inside(newpopn)
                            ## following line needed to return initial escapees 2021-11-23
                            newpopn[dispersing,] <- oldposition[dispersing,]
                            tries <- 0; maxtries <- 1000
                            while (any(dispersing) && tries <= maxtries) {
                                newpopn <- disperse(newpopn, turnoverpar, t,
                                    core, dispersing)  
                                # order of following rows reversed 2021-02-13
                                dispersing <- !inside(newpopn)
                                newpopn[dispersing,] <- oldposition[dispersing,]
                                tries <- tries+1
                            }
                        }
                        else if (turnoverpar$edgemethod == "none") {
                            ## no action, allow outward drift
                        }
                        else stop ("edgemethod ", turnoverpar$edgemethod, " not recognised")
                    }
                }
            }
            if (age) {
                covariates(newpopn)$age <- covariates(newpopn)$age + 1
            }
            gam <- turnoverpar$lambda[t] - turnoverpar$phi[t]
            if ((turnoverpar$recrmodel != 'specified') && (gam<0))
                stop ("invalid gamma in turnover")
            nrecruit <- switch (turnoverpar$recrmodel,
                constantN = nrow(oldpopn) - nsurv,
                discrete = discrete(gam * nrow(oldpopn)),
                binomial = rbinom(1, nrow(oldpopn), gam),
                poisson = rpois (1, gam * nrow(oldpopn)),
                multinomial = Nrecruits[t+1],
                specified = turnoverpar$Nrecruits[t],   # 2021-04-08
                -1)
            if (nrecruit<0) stop ("unrecognised recruitment model: ",turnoverpar$recrmodel)
            if (nrecruit>0) {
                
                recruits <- sim.popn(D = D, core = core, buffer = buffer,
                    model2D = model2D, buffertype = buffertype, poly = poly,
                    covariates = covariates, number.from = lastnumber + 1, 
                    Ndist = 'specified', Nbuffer = nrecruit,
                    nsessions = 1, details = details)
                # 2021-04-09
                if (age) {
                    if (is.null(covariates(recruits))) {
                        covariates(recruits) <- data.frame(age = rep(0, nrow(recruits))) 
                    }
                    else {
                        covariates(recruits)$age <- rep(0, nrow(recruits))
                    }
                }
                lastnumber <<- lastnumber + nrecruit
                newpopn <- rbind(newpopn, recruits, renumber = FALSE)
            }
            class(newpopn) <- class(MSpopn[[1]])
            attr(newpopn, 'mask') <- attr(MSpopn[[1]], 'mask')
            ## minor attributes are neglected... beware in future
            attr(newpopn, 'losses') <- nrow(oldpopn)-nsurv
            attr(newpopn, 'recruits') <- nrecruit
            newpopn
        }   ## end of turnover fn
        expands <- function (param, s) {
            if (is.null(param)) NULL
            else rep(param, length.out = s)
        }
        getbeta <- function (phi, lambda) {
            J <- length(phi)                 ## nsessions
            Nj <- cumprod(c(1,lambda[-J]))   ## relative number
            B <- Nj - c(0, (phi * Nj)[-J])   ## relative recruits
            B / sum(B)                       ## beta
        }
        
        #---------------------------------------------------------------------------------
        if (is.null(details$lambda)) {
            ## independent populations
            if (missing(D))
                D <- rep(NA, nsessions)
            else if (length(D) == 1)
                D <- rep(D, nsessions)
            else
                if (length(D) != nsessions) stop ("length(D) should equal nsessions")
            if (is.null(Nbuffer))
                Nbuffer <- rep(NA, nsessions)
            else if (length(Nbuffer) == 1)
                Nbuffer <- rep(Nbuffer, nsessions)
            else
                if (length(Nbuffer) != nsessions) stop ("length(Nbuffer) should equal nsessions")
            MSpopn <- mapply (session.popn, 1:nsessions, D, Nbuffer, Ndist, SIMPLIFY = FALSE)
        }
        else {
            ## projected population
            turnoverpar <- list(
                lambda = NULL, phi = 0.7, movemodel = 'static', sigma.m = 0, 
                move.a = NULL, move.b = 1, d0 = 0, edgemethod = "wrap",  
                survmodel = 'binomial', recrmodel = 'poisson', Nrecruits = 0)
            if (!is.null(details$wrap)) {
                warning("details option 'wrap' is deprecated; using edgemethod = 'wrap' if TRUE")
                if (details$wrap) details$edgemethod = 'wrap'
            }
            turnoverpar <- replace (turnoverpar, names(details), details)

            if (is.null(turnoverpar$move.a)) {
                ## allow legacy input as sigma.m
                turnoverpar$move.a <- turnoverpar$sigma.m
            }
            if (!is.function(turnoverpar$movemodel) &&
                    turnoverpar$movemodel == 'static' && turnoverpar$move.a != 0) {
                ## restore default from version < 3.2.1
                turnoverpar$movemodel <- 'BVN'
            }
            
            ################################################### 
            # 2021-10-15   
            # strip "zi", record for later
            if (is.function(turnoverpar$movemodel)) {
                turnoverpar$zeroinflated <- FALSE
            }
            else {
                turnoverpar$zeroinflated <- grepl('zi', turnoverpar$movemodel)
            }
            if (turnoverpar$zeroinflated) {
                turnoverpar$movemodel <- gsub("zi","",turnoverpar$movemodel)
                if (!(turnoverpar$movemodel %in% 
                        c('IND','BVN','BVE','RDE','uncorrelated','frE')))
                    stop ('zero-inflation not allowed with ', turnoverpar$movemodel)
            }
            ################################################### 
            
            turnoverpar$lambda  <- expands(turnoverpar$lambda, nsessions)
            turnoverpar$phi     <- expands(turnoverpar$phi, nsessions)
            turnoverpar$move.a <- expands(turnoverpar$move.a, nsessions)
            turnoverpar$move.b <- expands(turnoverpar$move.b, nsessions)
            turnoverpar$Nrecruits <- expands(turnoverpar$Nrecruits, nsessions-1)
            MSpopn <- vector(nsessions, mode = 'list')
         
            if (turnoverpar$recrmodel == "multinomial") {
                beta <- getbeta(turnoverpar$phi, turnoverpar$lambda)
                if (is.null(details$superN)) {
                    ## infer superN from a trial initial population
                    ## this may introduce some variation
                    N <- nrow(session.popn(1, D, Nbuffer, Ndist))
                    details$superN <- N/beta[1]
                    warning("multinomial recruitment superN not specified, using ", 
                            round(details$superN,1))
                }
                Nrecruits <- rmultinom(1, details$superN, beta)
                MSpopn[[1]] <- session.popn(1, D, Nrecruits[1], Ndist="specified")
            }
            else {
                MSpopn[[1]] <- session.popn(1, D, Nbuffer, Ndist)
            }
            lastnumber <<- lastnumber + nrow(MSpopn[[1]])
            # 2021-04-09
            if (age) {
                if (is.null(covariates(MSpopn[[1]])))
                    covariates(MSpopn[[1]]) <- data.frame(age = rep(0,nrow(MSpopn[[1]])))
                else 
                    covariates(MSpopn[[1]])$age <- rep(0,nrow(MSpopn[[1]]))
            }
            for (i in 2:nsessions) {
                MSpopn[[i]] <- turnover(MSpopn[[i-1]], i-1)
            }
        }
        if (model2D == 'linear')
            class(MSpopn) <- c('linearpopn', 'popn', 'list')
        else
            class(MSpopn) <- c('popn', 'list')
        names(MSpopn) <- 1:nsessions
        MSpopn
    }
    else {
        ##########################
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
        ##########################

        # 2022-06-06 simplify
        # getnm <- function (scaleattribute = 'area', scale = 1, D) {
        getnm <- function (cellsize, D) {
            ## 2014-09-03 for "IHP" and "linear"
            if ((length(D) == 1) & (is.character(D))) {
                D <- covariates(core)[,D]
            }
            else if ((length(D) == 1) & (is.numeric(D))) {
                D <- rep(D, nrow(core))
            }
            
            if (any(is.na(D))) {
                D[is.na(D)] <- 0
                warning ("NA values of D set to zero")
            }
            if (any(D<0)) {
                D <- pmax(0,D)
                warning ("negative D set to zero")
            }
            ## D vector, 1 per cell
            D <- D * cellsize
            if (!is.null(Nbuffer)) {  ## includes Ndist == 'specified'
              N <- round(Nbuffer)
            }
            else {
              N <- sum(D)
            }
            if (Ndist == 'poisson') {
                N <- rpois(1,N)
            }
            # 2023-05-30 catch boundary case
            if (N<=0) {   
                rep(0,length(D))
            }
            else {
                rmultinom (1, N, D)
            }
        }
        ##########################

        ## 2016-03-08
        ## allow core to be a named object
        if (is.character(core)) core <- get(core)
            
        if (model2D %in% c('IHP')) {
            if (!inherits(core, 'mask'))
                stop ("for model2D = IHP, 'core' should be a habitat mask")
            # typo fixed 2023-09-17
            if (nsessions>1 && details$movemodel!='static' && 
                    ('settle' %in% names(covariates(core))) && 
                    details$edgemethod %in% c('truncate', 'normalize')) 
            {
                s <- covariates(core)$settle
                if (any(is.na(s)) || min(s)<0 || max(s)>1) 
                    stop('settle covariate should be in range 0-1 with none missing')
            }
            cellsize <- attr(core, 'area')
            
            ## 2022-11-24 function D
            if (is.function(D)) {
                D <- D(mask = core, parm = details)
                covariates(core)$D <- D 
            }

            nm <- getnm(cellsize, unlist(D))
            jitter <- matrix ((runif(2*sum(nm))-0.5) * attr(core,'spacing'), ncol = 2)
            animals <- core[rep(1:nrow(core), nm),] + jitter
            animals <- as.data.frame(animals)
            xl <- range(animals[,1])
            yl <- range(animals[,2])
        }
        else  if (model2D == 'linear') {
            if (!inherits(core, 'linearmask'))
                stop ("for model2D = linear, 'core' should be a linear mask")
            cellsize <- attr(core, 'spacing') * 0.001
            nm <- getnm(cellsize, D)
            animals <- core[rep(1:nrow(core), nm),] * 1  ## * 1 to shed attributes...
            animals <- as.data.frame(animals)
            xl <- range(animals[,1])
            yl <- range(animals[,2])
            
        }
        else {
            # population in arena +/- buffer from traps
            # 2023-11-06 buffer should equal half spacing if (inherits(core, 'mask'))
            buff <- c(-buffer,+buffer)
            xl <- range(core$x) + buff
            yl <- range(core$y) + buff
            area <- diff(xl) * diff(yl) * 0.0001  # ha not sq metres

            bufferpoly <- switch(buffertype,
                                 rect = NA,
                                 convex = buffer.contour(core, buffer = buffer,
                                                         convex = TRUE, plt = FALSE)[[1]],
                                 concave = buffer.contour(core, buffer = buffer,
                                                          convex = FALSE, plt = FALSE))

            bufferarea <- switch (buffertype,
                                  rect = area,
                                  convex = polyarea (bufferpoly),
                                  concave = sum(sapply(bufferpoly, polyarea)))
            
            if ((buffertype == 'concave') & (is.null(Nbuffer)))
                warning("automatic Nbuffer unreliable with concave buffer")

            ## If target number (Nbuffer) not specified, get from density x area
            ## N is the number of centres to simulate in the unbuffered (rectangular) area
            if (is.null(Nbuffer)) {
                N  <- switch (Ndist,
                              poisson = rpois(1, lambda = D[1] * area),
                              fixed = discreteN (1, D[1] * area),
                              specified = round(D[1]))
                Nbuffer  <- switch (Ndist,
                                    poisson = rpois(1, lambda = D[1] * bufferarea),
                                    fixed = discreteN (1, D[1] * bufferarea),
                                    specified = round(D[1]))
            }
            else {
                Nbuffer  <- switch (Ndist,
                                    poisson = rpois(1, lambda = Nbuffer),
                                    fixed = discreteN (1, Nbuffer),
                                    specified = round(Nbuffer))
                N <- Nbuffer
            }

            if (model2D == 'poisson') {
                animals <- data.frame (
                    x = runif(N)*diff(xl)+xl[1],
                    y = runif(N)*diff(yl)+yl[1])

                ## 2014-12-29, 2015-01-11, 2015-01-13 allow buffering / poly
                if (buffertype %in% c('convex','concave')) {

#                     if (Ndist == 'fixed') {
                        maxtries <- 10
                        tries <- 1
                        animals <- switch(buffertype,
                                      convex = animals[pointsInPolygon(animals, bufferpoly),],
                                      concave = animals[distancetotrap(animals, core)<= buffer,])
                        ## repeat if not enough
                        while ((nrow(animals) < Nbuffer) & (tries < maxtries)) {
                            animals <- rbind(animals, data.frame (
                                                                  x = runif(N)*diff(xl)+xl[1],
                                                                  y = runif(N)*diff(yl)+yl[1]))
                            animals <- switch(buffertype,
                                        convex = animals[pointsInPolygon(animals, bufferpoly),],
                                        concave = animals[distancetotrap(animals, core)<= buffer,])
                            tries <- tries + 1
                        }
                        if (tries >= maxtries)
                            warning("exceeded maxtries in sim.popn")
                        animals <- animals[1:Nbuffer,]
#                     }
#                     else {
#                         animals <- switch(buffertype,
#                                       convex = animals[pointsInPolygon(animals, bufferpoly),],
#                                       concave = animals[distancetotrap(animals, core)<= buffer,])
#                     }
                }
            }
            else if (model2D == 'coastal') {
                if (is.null(details$Beta))
                    details$Beta <- c(1,1.5,5,1)
                a1 <- details$Beta[1]
                b1 <- details$Beta[2]
                a2 <- details$Beta[3]
                b2 <- details$Beta[4]
                animals <- data.frame (x = rbeta(N,a1,b1)*diff(xl)+xl[1],
                                       y = rbeta(N,a2,b2)*diff(yl)+yl[1])
            }
            else if (model2D=='hills') {
                hills <- details$hills
                if (is.null(hills)) hills <- c(1,1)
                hills <- c(hills, rep(0,4-length(hills)))
                nhillx <- abs(hills[1])
                nhilly <- abs(hills[2])
                offset <- any(hills[1:2]<0)
                dx <- hills[3]
                dy <- hills[4]
                dx <- ifelse (dx<0, runif(1), dx)
                dy <- ifelse (dy<0, runif(1), dy)
                xhill <- sample(0:(nhillx-1), N, replace=TRUE)
                yhill <- sample(0:(nhilly-1), N, replace=TRUE)
                x <- asin(runif(N)*2 - 1)/pi + 0.5 + xhill
                y <- asin(runif(N)*2 - 1)/pi + 0.5 + yhill
                if (offset) x <- x + (yhill %% 2) * 0.5
                x <- x/nhillx + dx
                y <- y/nhilly + dy
                x <- ifelse (x>1, x-1, x)
                y <- ifelse (y>1, y-1, y)
                animals <- data.frame (x = x * diff(xl)+xl[1],
                                       y = y * diff(yl)+yl[1])
            }
            else if (model2D == 'cluster') {
                ## Neyman-Scott distribution with wrapping
                xrange <- diff(xl)
                yrange <- diff(yl)
                if (details$mu<=0) {
                    nparent <- N   ## not clustered
                    offspr <- sweep(matrix(runif(2*nparent), ncol = 2), 2, c(xrange,yrange), '*')
                }
                else {
                    nparent <- switch (Ndist,
                                       poisson = rpois(1, lambda=D[1] * area/details$mu),
                                       fixed = discreteN (1, D[1] * area / details$mu),
                                       specified = discreteN (1, D[1] / details$mu))  ## here arg D is N
                    if (nparent==0)
                        warning ("zero clusters")
                    parent <-  sweep(matrix(runif(2*nparent), ncol = 2), 2, c(xrange,yrange), '*')
                    # number of offspring for each parent
                    if (!is.null(details$clone) && details$clone == 'constant') {
                        noffspr <- rep(details$mu, nparent)
                    }
                    else {
                        noffspr <- rpois(nparent, details$mu)
                    }
                    N <- sum(noffspr)
                        # for backward compatibility
                    if (is.null(details$scale)) details$scale <- details$hsigma
                    # scale = 0 to clone parent locations with no displacement 2023-10-10
                    offspr <- matrix(rnorm(2*N), ncol = 2) * details$scale
                    if (N>0) {
                        parentn <- rep(1:nparent, noffspr)
                        offspr <- offspr + parent[parentn,,drop = FALSE]
                        # toroidal wrapping
                        while (any ((offspr[,1]<0) | (offspr[,1]>xrange) | (offspr[,2]<0) |
                                    (offspr[,2]>yrange))) {
                            offspr[,1] <- ifelse (offspr[,1]<0, offspr[,1]+xrange, offspr[,1])
                            offspr[,1] <- ifelse (offspr[,1]>xrange, offspr[,1]-xrange, offspr[,1])
                            offspr[,2] <- ifelse (offspr[,2]<0, offspr[,2]+yrange, offspr[,2])
                            offspr[,2] <- ifelse (offspr[,2]>yrange, offspr[,2]-yrange, offspr[,2])
                        }
                    }
                }
                animals <- as.data.frame(sweep(offspr,2,c(xl[1],yl[1]),'+'))
            }
            else if (model2D == 'even') {
                ## 'even' distribution from Efford 2004 and Density
                D <- N / area
                xrange <- diff(xl)
                yrange <- diff(yl)
                cellside <- sqrt(10000/D)
                centrex <- seq(xl[1]+cellside/2, xl[2]+3*cellside/2, cellside)
                centrey <- seq(yl[1]+cellside/2, yl[2]+3*cellside/2, cellside)
                centres <- expand.grid(x=centrex, y=centrey)
                animals <- centres + runif(nrow(centres*2)) * cellside - cellside/2
                animals <- animals[animals[,1] <= xl[2] & animals[,2] <= yl[2], ]
                # possibly save grid?
            }
            else  if (model2D %in% c("rLGCP", "rThomas")) {
                if (requireNamespace("spatstat.geom", quietly = TRUE) && 
                    requireNamespace("spatstat.random", quietly = TRUE)) {
                    if (!is.numeric(D) || length(D)>1) {
                        stop ("for model2D in (rLGCP, rThomas) D should be a scalar")
                    }
                    if (is.null(details$saveLambda)) details$saveLambda <- TRUE
                    # spatstat window
                    ow <- spatstat.geom::owin(xl, yl)  
                    if (is.null(details$eps)) details$eps <- diff(xl)/64
                    if (model2D == 'rLGCP') {
                        # D, var, scale
                        # mu for rLGCP is derived from D, var
                        mu <- log(D/1e4) - details$var/2    # mean density / m^2 on log scale
                        pts <- spatstat.random::rLGCP(
                            model      = "exp",
                            mu         = mu,
                            var        = details$var,
                            scale      = details$scale,
                            win        = ow,
                            saveLambda = details$saveLambda,
                            eps        = details$eps)
                    }
                    else if (model2D == 'rThomas') {
                        # kappa for rThomas is D/mu
                        kappa <- D/1e4/details$mu   # mean density / m^2
                        pts <- spatstat.random::rThomas(
                            kappa       = kappa,
                            scale       = details$scale,
                            mu          = details$mu,
                            win         = ow,
                            nonempty    = FALSE, 
                            saveparents = TRUE,
                            saveLambda  = details$saveLambda,
                            eps         = details$eps)    # 2023-11-07
                    }
                    animals <- spatstat.geom::coords(pts)
                    animals <- as.data.frame(animals)
                    if (details$saveLambda) {
                        attr(animals, "Lambda") <- im2mask(attr(pts, "Lambda")) 
                    }
                    if (model2D == 'rThomas') {
                        attr(animals, "parents") <- as.data.frame(attr(pts, "parents"))
                    }
                }
                else {
                    stop ("rLGCP and rThomas use the package spatstat")
                }
            }
            else stop ("unrecognised 2-D distribution")
        }
        names(animals) <- c('x','y')
        attr(animals,'covariates') <- NULL
        if (!is.null(covariates)) {
            tempcov <- list()
            for (i in 1:length(covariates)) {
               covi <- sample (names(covariates[[i]]), replace = T, size= nrow(animals),
                               prob=covariates[[i]])
               temptxt <- paste ('tempcov$', names(covariates[i]), '<- covi',
                               sep = '')
               eval(parse(text=temptxt))
            }
            attr(animals,'covariates') <- as.data.frame(tempcov)
        }
        if (nrow(animals) > 0)   ## condition added 2011-03-27
            row.names (animals) <- number.from : (nrow(animals)+number.from-1)

        if (keep.mask && model2D %in% c('IHP','linear')) {
                attr(animals, 'mask') <- core
        }
        if (model2D == 'linear') {
            class(animals) <- c('linearpopn', 'popn', 'data.frame')
        }
        else {
            class(animals) <- c('popn', 'data.frame')
            ##-------------------------
            ## restrict to a polygon
            ## added 2011-10-20
            ## NOTE 2014-12-29 this breaks Ndist = 'fixed'
            if (!is.null(poly)) {
                animals <- subset(animals, poly = poly, ...)
            }
            ##-------------------------
        }
        attr(animals, 'seed') <- RNGstate   ## save random seed
        attr(animals, 'Ndist') <- Ndist
        attr(animals, 'Nbuffer') <- Nbuffer
        attr(animals, 'D') <- D   # bulky if mask
        attr(animals, 'model2D') <- model2D
        attr(animals, 'buffertype') <- buffertype
        attr(animals, 'boundingbox') <- expand.grid (x=xl,y=yl)[c(1,3,4,2),]
        animals
    }
}
