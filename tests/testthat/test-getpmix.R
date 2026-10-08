# mixture proportions by animal: the first detector-occasion in use need not be
# detector 1 on occasion 1 (e.g. detectors installed late); secr 5.5.1

library(secr)

## PIA: 1 x animals x occasions x detectors x classes. Only detectors 3-4 on occasions 2-3
## are in use, so the first used (s,k) is occasion 2 of detector 3 (flattened index 8 of 12)
S <- 3; K <- 4; nc <- 2; nmix <- 2
PIA <- array(0L, c(1, nc, S, K, nmix))
PIA[1, , 2:3, 3:4, 1] <- 1L          # class 1 -> row 1 of realparval
PIA[1, , 2:3, 3:4, 2] <- 2L          # class 2 -> row 2 of realparval
realparval <- cbind(pmix = c(0.3, 0.7))

test_that("secr_firstsk indexes occasion fastest", {
    expect_equal(as.vector(secr:::secr_firstsk(PIA[1, 1, , , 1, drop = FALSE])), 8)
})

test_that("secr_getpmix with detector 1 unused on early occasions", {
    ## animal 1 is known to be in class 1 (knownclass 2); animal 2 has unknown class (1)
    pmixn <- secr:::secr_getpmix(c(2, 1), PIA, realparval)
    expect_equal(as.vector(pmixn[, 1]), c(1, 0))
    expect_equal(as.vector(pmixn[, 2]), c(0.3, 0.7))
    expect_equal(as.vector(attr(pmixn, "pmix")), c(0.3, 0.7))
})

test_that("getpmixall (sim.detect) with detector 1 unused on early occasions", {
    expect_equal(secr:::getpmixall(PIA, realparval), c(1, 2))
})
