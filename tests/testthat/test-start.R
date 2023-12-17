# test makeStart

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
system.time({
msk2 <- make.mask(traps(captdata), buffer = 100, spacing = 20)
set.seed(123)
fit0 <- secr.fit(captdata, 
                 mask  = msk2, 
                 model = sigma~b, 
                 trace = FALSE,
                 ncores = 2,
                 details = list(fixedbeta = c(NA,NA,NA,0.31)))
parindx0 <- list(g0=1, sigma=2)
parindx1 <- list(D=1, g0=2, sigma=3)
parindx2 <- list(D=1, g0=2, sigma=c(3,4))
links <- list(D='log', g0='logit', sigma='log')

test_that("correct start vector", {
    
    S0 <- makeStart(start    = secrdemo.CL, 
                    parindx  = parindx0)
    
    S1 <- makeStart(start    = secrdemo.0, 
                    parindx  = parindx2,
                    details  = list(fixedbeta = c(NA,NA,NA,0.31)))
    
    S2 <- makeStart(start    = fit0, 
                    parindx  = fit0$parindx,
                    details  = list(fixedbeta = c(NA,NA,NA,0.31)))
    
    S3 <- makeStart(start    = fit0, 
                    parindx  = fit0$parindx,
                    details  = list(fixedbeta = NULL))
    
    set.seed(123) # for autoini
    S4 <- makeStart(start    = NULL, 
                    parindx  = parindx1,
                    capthist = captdata, 
                    mask     = msk2,
                    detectfn = 0,
                    link     = links,
                    details  = list(fixedbeta = NULL, trace = FALSE))
    
    expect_equal(S0, c(-0.9784944, 3.3798320), tolerance = 1e-4, 
                 check.attributes = FALSE)
    expect_equal(S1, c(1.7010692, -0.9784937, 3.3798317), tolerance = 1e-4, 
                 check.attributes = FALSE)
    expect_equal(S2, c(1.873340, -1.144755, 3.154245), tolerance = 1e-4, 
                 check.attributes = FALSE)
    expect_equal(S3, c(1.873340, -1.144755, 3.154245, 0.3100), tolerance = 1e-4, 
                 check.attributes = FALSE)
    expect_equal(S4, c(1.717297, -1.182687, 3.426041), tolerance = 1e-4, 
                 check.attributes = FALSE)
})

})
