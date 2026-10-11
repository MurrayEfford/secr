# log-likelihood must not depend on details$uselog, including when sigma is small relative
# to the mask so that the hazard is effectively zero in many mask cells; secr 5.5.1
# (a cell with negligible hazard on an occasion without capture was wrongly given probability 0)

library(secr)

test_that("uselog = TRUE and FALSE give the same log-likelihood (small sigma, wide mask)", {
    ch <- ovenCH[[1]]
    msk <- make.mask(traps(ch), buffer = 300, nx = 24, type = "trapbuffer")
    LL <- function(CL, uselog) {
        secr.fit(ch, mask = msk, CL = CL, detectfn = "HN", trace = FALSE,
                 start = if (CL) list(g0 = 0.3, sigma = 30) else list(D = 1, g0 = 0.3, sigma = 30),
                 details = list(LLonly = TRUE, uselog = uselog))
    }
    for (CL in c(TRUE, FALSE)) {
        expect_equal(as.numeric(LL(CL, FALSE)), as.numeric(LL(CL, TRUE)), tolerance = 1e-8)
    }
    ## reference values from secr 5.4.3 (no uselog option)
    expect_equal(as.numeric(LL(FALSE, FALSE)), -219.7162, tolerance = 1e-6)
    expect_equal(as.numeric(LL(TRUE, FALSE)), -216.6559, tolerance = 1e-6)
})
