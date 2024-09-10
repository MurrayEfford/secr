# test MCgof

library(secr)

test_that("MCgof correct with simple multicatch model", {
    # switch to multi from single generates warning
    expect_warning(test0 <- MCgof(secrdemo.0, nsim=2, seed=123, quiet = TRUE))
    discrepancytable <- summary(test0)
    target1 <- c(yik = 468.87675, yi = 207.48601, yk = 137.77262)
    target2 <- c(yik = 459.06023, yi = 205.23678, yk = 150.62346)
    # Tobs, Tsim rows of summary
    expect_equal(discrepancytable['Tobs',], target1) 
    expect_equal(discrepancytable['Tsim',], target2) 
})
