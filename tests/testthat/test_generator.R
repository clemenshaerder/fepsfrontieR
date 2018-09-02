context ("sfmfep")

library (MASS)
library (dplyr)
library (truncnorm)
library (dplyr)
library (stats)
library (numDeriv)
library (parallel)
library (magrittr)
library (stringr)

test_that ("sfmfep works", {

  # Test input
  N <- 40;
  Time <- 5
  beta <- c(0.5)
  delta <- c(0.5)
  sigma_u <- 0.2
  sigma_v <- 0.1
  t.formula <- formula(y  ~ x + (z))
  test.data <- sfm.data  # package data
  method <- "within"
  boot <- F
  mu = 0
  myPar = NULL
  boot = F
  B = NULL
  alphaCI <- 0.05
  estimate = T
  panel = NULL
  parallel = F

  # tests if bootstrapping works for method = "firstdiff"
  summary(firstdiffBoot)
  firstdiffBoot <- sfmfep(formula = t.formula, bootstrap = T, B = 30, alphaCI = c(0.1, 0.05),
                          method = method,  N = N, Time = Time, parallel = parallel,
                          data = test.data, mu = mu, myPar = myPar)
  expect_type (object = firstdiffBoot, type = "list")

  # tests if bootstrapping works for method = "within"
  withinBoot <- sfmfep(formula = t.formula, bootstrap = F, B = 10, parallel = parallel,
                       method = "firstdiff", N = N, Time = Time, data = test.data,
                       mu = mu, myPar = myPar)
  expect_type (object = withinBoot, type = "list")

  # Tests if optim " N & T" works & mu > 0
  balancedNT <- sfmfep(formula = t.formula, N = 40, Time = 5, bootstrap = boot, B,
                       method = method, parallel = parallel,
                       data = test.data, mu = 1, myPar = myPar)
  expect_type (object = balancedNT, type = "list")

  # unbalanced & incorrect data dimension test
  expect_error ( sfmfep (formula = t.formula, method = method,
                         bootstrap = boot, B = B, parallel = parallel,
                         N=40,Time=c(10,9,8,3,1,6,7,7,8,9), data = test.data, mu = mu, myPar = myPar))

  # Tests if option "panel" works when we specify the column.
  panelColTest1 <- sfmfep (formula = t.formula, method = method,
                           bootstrap = boot, B = B, parallel = parallel,
                           panel = test.data$producer, data = test.data, mu = mu, myPar = myPar)
  expect_type (object = panelColTest1, type = "list")

  # Tests if option "panel" works when we specify the column.
  panelColTest2 <- sfmfep (formula = t.formula, method = method,
                           bootstrap = boot, B = B, parallel = parallel,
                           panel = test.data[, 1], data = test.data, mu = mu, myPar = myPar)
  expect_type (object = panelColTest2, type = "list")

  # Tests if option "panel" works
  panelTest <- sfmfep (formula = t.formula, method = method,
                       bootstrap = boot, B = B, parallel = parallel,
                       panel = "producer", data = test.data, mu = mu, myPar = myPar)
  expect_type (object = panelTest, type = "list")

  # Tests if defined starting points "myPar" works with Bootstrapping & panel
  panelMyParBoot <- sfmfep(formula = t.formula, method = method, panel ="producer", bootstrap = T, B = 5,
                           data = test.data, mu = mu, parallel = parallel,
                           myPar = c(sigma_u = 1, sigma_v=2, beta = c(1), delta = c(1)))
  expect_type (object = panelMyParBoot, type = "list")

  # Tests if it works when estimates are provided (estimate = F)
  panelEstimateF <- sfmfep(formula = t.formula, method = method, panel ="producer", bootstrap = boot, B = B,
                       data = test.data, mu = mu, estimate = F, parallel = parallel,
                       myPar = c(sigma_u = 1, sigma_v = 2, beta = c(1), delta = c(2)))
  expect_type (object = panelEstimateF, type = "list")

  # Tests if it works when CIs are not wanted
  panelNoCI <- sfmfep(formula = t.formula, method = method, panel = "producer", bootstrap = boot, B = B,
                       data = test.data, mu = mu, alphaCI = NULL, parallel = parallel,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(2), delta = c(1)))
  expect_type (object = panelNoCI, type = "list")

  # Tests unbalanced panels with panels without CI
  test.data <- test.data[-60, ]
  unbalancedpanel <- sfmfep(formula = t.formula, method = method, panel ="producer", bootstrap = boot, B = B,
                       data = test.data, mu = mu, alphaCI = 0.05, parallel = parallel,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(2), delta = c(1)))
  expect_type (object = unbalancedpanel, type = "list")

  unbalancedBoot <- sfmfep(formula = t.formula, method = method, panel ="producer", bootstrap = T, B = 5,
                            data = test.data, mu = mu, alphaCI = NULL, parallel = parallel,
                            myPar = c(sigma_u = 1, sigma_v=2, beta = c(2), delta = c(1)))
  expect_type (object = unbalancedBoot, type = "list")

  unbalancedBootNT <- sfmfep(formula = t.formula, method = method, N = 40, Time =c(rep(5,39),4), bootstrap = T, B = 5,
                           data = test.data, mu = mu, alphaCI = 0.05, parallel = parallel,
                           myPar = c(sigma_u = 1, sigma_v=2, beta = c(2), delta = c(1)))
  expect_type (object = unbalancedBoot, type = "list")

  rm(list=ls(all=TRUE))
})


context ("SFM.generate")
library (truncnorm)

test_that ("SFM.generate creates a correct output format", {

  N <- 2
  Time <- 5
  beta <- c(0.5, 0.5)
  delta <- c(0.5, 0.5)
  sigma_u = 0.1
  sigma_v = 0.2
  mu <- 0.5

  testSFM1 <- SFM.generate(N = N, Time = Time, beta = beta, delta = delta,
                           mu = mu, sigma_u = sigma_u, sigma_v = sigma_v)
  expect_type (object = testSFM1, type = "list")

  N <- "c(3,2)"
  expect_error( SFM.generate(N = N, Time = Time, beta = beta, delta = delta,
                             mu = mu, sigma_u = sigma_u, sigma_v = sigma_v))

  rm(list=ls(all=TRUE))
})
