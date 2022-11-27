#

supplement.grid <- function (traps, spfactor = 0.2, randomfraction = NULL, skip = 0, plt = FALSE, ...) {
    if (ms(traps)) stop ("traps should be single-session")
    if (spfactor>= 0.5) stop ("spfactor should be less than 0.5")
    sp <- spacing(traps)
    if (is.null(randomfraction)) {
        skip <- skip+1
        nx <- round((length(unique(traps$x))-1)/skip) + 1
        ny <- round((length(unique(traps$y))-1)/skip) + 1
        newpoints <- make.grid(nx = nx, ny = ny, spacing = sp*skip, 
            originxy = apply(traps,2,min), detector = detector(traps))
    }
    else {
        newpoints <- subset(traps, sample.int(nrow(traps), randomfraction * nrow(traps)))
    }
    theta <- seq(0,3*pi/2,pi/2)[sample.int(4,nrow(newpoints), replace = TRUE)]
    dx <- cos(theta)*spfactor*sp
    dy <- sin(theta)*spfactor*sp
    newpoints[,] <- newpoints + cbind(dx, dy)
    if (plt) {
        plot(traps, ...)
        dots <- list(...)
        dp <- dots[['detpar']]
        dp$pch <- 16
        args <- c(list(x = newpoints, add = TRUE, detpar = dp), dots[names(dots) != 'detpar'])
        do.call(plot, args)
    }
    rbind(traps, newpoints)
}
# 
# par(mfrow=c(1,2))
# 
# g1 <- make.grid(7,7, spacing = 50)
# g2 <- supplement.grid(g1, skip = 2, plt = TRUE, border=30, gridlines = FALSE)
# 
# g3 <- supplement.grid(g1, randomfraction=0.5, plt=T, border=30, gridlines = FALSE)
# g4 <- make.grid(7,11, spacing = 50)
# 
# 
# library(secrdesign)
# scen <- make.scenarios(D=10, lambda0=0.2, noccasions=5, sigma = 20, trapsindex=1:3, detectfn = 'HHN')
# scenarioSummary(scen, trapset=list(g1,g3,g4))

# $1
# 1 
# n     mean      se
# estimate    10  9.62526 0.30719
# SE.estimate 10  1.53513 0.08917
# lcl         10  7.05871 0.19367
# ucl         10 13.13686 0.50988
# RB          10 -0.03747 0.03072
# RSE         10  0.15868 0.00509
# COV         10  1.00000 0.00000
# 
# $2
# 2 
# n     mean      se
# estimate    10  9.92822 0.39823
# SE.estimate 10  1.17573 0.04281
# lcl         10  7.87845 0.32793
# ucl         10 12.51237 0.48601
# RB          10 -0.00718 0.03982
# RSE         10  0.11866 0.00151
# COV         10  1.00000 0.00000
# 
# $3
# 3 
# n     mean      se
# estimate    10  9.99275 0.27035
# SE.estimate 10  1.22247 0.04665
# lcl         10  7.87016 0.20423
# ucl         10 12.69019 0.36659
# RB          10 -0.00073 0.02704
# RSE         10  0.12215 0.00234
# COV         10  1.00000 0.00000
