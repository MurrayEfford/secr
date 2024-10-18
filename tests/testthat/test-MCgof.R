# test MCgof, summary.MCgof

library(secr)

# test gives inconsistent results on different platforms
# test_that("MCgof correct with simple multicatch model", {
#     # switch to multi from single generates warning
#     expect_warning(test0 <- MCgof(secrdemo.0, nsim=2, seed = 123, quiet = TRUE))
#     discrepancytable <- summary(test0)
#     target1 <- c(yik = 513.28664, yi = 252.48601, yk = 158.46254)
#     target2 <- c(yik = 504.04813, yi = 240.39090, yk = 165.88048)
#     # Tobs, Tsim rows of summary
#     expect_equal(discrepancytable['Tobs',], target1, tolerance = 0.001) 
#     expect_equal(discrepancytable['Tsim',], target2, tolerance = 0.001) 
# })

test_that("MCgof correct with proximity model", {
    expect_warning(test0 <- MCgof(blackbear.0, nsim=2, seed = 123, quiet = TRUE))
    discrepancytable <- summary(test0)
    target1 <- c(yik = 316.43949, yi = 99.12776, yk = 17.80118)
    target2 <- c(yik = 291.53224, yi = 85.18867, yk = 21.97622)
    # Tobs, Tsim rows of summary
    expect_equal(discrepancytable['Tobs',], target1, tolerance = 0.001) 
    expect_equal(discrepancytable['Tsim',], target2, tolerance = 0.001) 
})
