## 2025-08-03
## tests of group factor based on multiple individual covariates
## (the 'groups' argument of secr.fit, derived.secr)

library(secr)

grouplevels <- secr:::secr_group.levels(OVpossumCH[[1]], c('Sex','Ageclass'))
factorlevels <- levels(secr:::secr_group.factor(OVpossumCH[[1]], c('Sex','Ageclass')))

# detection models
detectiondesign <- secr.design.MS (OVpossumCH[[1]], model = list(lambda0=~g), 
                                   groups = c('Sex','Ageclass'))
detectiondesignlevels <- colnames(detectiondesign[[1]][[1]])

# density models
expect_warning(testCH <- subset(OVpossumCH[[1]], 180:223, trap=31:60))
expect_warning(fit <- secr.fit(testCH, CL = FALSE, model=list(D~g), detectfn = 'HEX', trace = FALSE,
                               start = c(-0.615, 0.575, 2.314, 2.476, -1.224, 2.787),
                               groups = c('Sex','Ageclass'), method = 'none'))
densitydesignlevels <- colnames(fit$designD)

# new data
newdatalevels <- levels(makeNewData(fit)$g)             

baselevels <- c("F.1", "M.1", "F.2", "M.2")
test_that("group levels OK", {
    expect_equal(grouplevels, baselevels)
    expect_equal(factorlevels, baselevels)
    expect_equal(detectiondesignlevels, c("lambda0.(Intercept)", "lambda0.gM.1", 
                                          "lambda0.gF.2", "lambda0.gM.2"))
    expect_equal(densitydesignlevels, c("(Intercept)", "gM.1", "gF.2", "gM.2"))
    expect_equal(newdatalevels, baselevels)
})

