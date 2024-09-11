# test MCgof, summary.MCgof

library(secr)

test_that("MCgof correct with simple multicatch model", {
    # switch to multi from single generates warning
    expect_warning(test0 <- MCgof(secrdemo.0, nsim=2, seed = 123, quiet = TRUE))
    discrepancytable <- summary(test0)
    target1 <- c(yik = 462.05233, yi = 199.98601, yk = 130.30719)
    target2 <- c(yik = 468.53801, yi = 209.66753, yk = 148.41228)
    # Tobs, Tsim rows of summary
    expect_equal(discrepancytable['Tobs',], target1, tolerance = 0.001) 
    expect_equal(discrepancytable['Tsim',], target2, tolerance = 0.001) 
})