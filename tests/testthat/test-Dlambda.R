# test Dlambda model
# 2026-01-20

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

msk <- make.mask(traps(ovenCH[[1]]), buffer = 300, type='trapbuffer', nx = 32)

test_that("Dlambda likelihood correct", {
    LL1 <- secr.fit(ovenCH, model = D~1, mask = msk, trace = TRUE, 
                    details = list(Dlambda = TRUE, LLonly = TRUE), 
                    start = c(0,0,-3.5,4), ncores = 2)
    expect_equal(LL1, -956.66706, tolerance = 1e-4, check.attributes = FALSE)
})

# possibly include Dlambda fitted model in ovenbird, and
#predictDlambda(fit1)
