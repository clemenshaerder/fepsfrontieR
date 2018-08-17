context ("sfmfep")

library (MASS)
library (dplyr)

test_that ("sfmfep works", {

  # Test input
  N <- 2;
  Time <- 30
  beta <- c(0.5,3);
  delta <- c(0.5,2)
  sigma_u <- 2;
  sigma_v <- 1
  t.formula <- formula(y  ~ x1 + x2 + (z1 + z2))
  test.data <- sfm.data  # package data
  method <- "within"
  boot <- F
  mu = 0
  myPar = NULL
  boot = F
  B = NULL


  # tests if bootstrapping works for method = "firstdiff"
  firstdiffBoot <- sfmfep(formula = t.formula, bootstrap = T, B = 10, method = method,
                          N = N, Time = Time, data = test.data, mu = mu, myPar = myPar)
  expect_type (object = firstdiffBoot, type = "list")

  # tests if bootstrapping works for method = "within"
  withinBoot <- sfmfep(formula = t.formula, bootstrap = T, B = 10, method = "within",
                       N = N, Time = Time, data = test.data, mu = mu, myPar = myPar)
  expect_type (object = withinBoot, type = "list")

  # Tests if optim " N & T" works & mu > 0
  balancedNT <- sfmfep(formula = t.formula, N = 2,Time = 30, bootstrap = boot, B, method = method,
                       data = test.data, mu = 1, myPar = myPar)
  expect_type (object = balancedNT, type = "list")

  # unbalanced & incorrect data dimension test
  expect_error ( sfmfep (formula = t.formula, method = method, bootstrap = boot, B = B,
                         N=2,Time=c(30,29), data = test.data, mu = mu, myPar = myPar) )

  # Tests if option "group" works TODO:(currently throws a NaN but result is correct)
  groupTest <- sfmfep (formula = t.formula, method = method, bootstrap = boot, B = B,
                       group = "gr", data = test.data, mu = mu, myPar = myPar)
  expect_type (object = groupTest, type = "list")

  # Tests if defined starting points "myPar" works with Bootstrapping & group
  GroupMyParBoot <- sfmfep(formula = t.formula, method = method, group ="gr", bootstrap = T, B = 5,
                       data = test.data, mu = mu,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
  expect_type (object = GroupMyParBoot, type = "list")

  # Tests if it works when estimates are provided (estimate = F)
  # TODO(Olli) -> This test doenst provide an output. I guess it is due to the negative Fisher of sigma_v
  # Anyhow, a CI is computed for the other parameters. Thus, there should still be an ouput.

  GroupEstimateF <- sfmfep(formula = t.formula, method = method, group ="gr", bootstrap = boot, B = B,
                       data = test.data, mu = mu, estimate = F,
                       myPar = c(sigma_u = 1, sigma_v = 2, beta = c(1,2), delta = c(1, 2)))
  expect_type(object = GroupEstimateF, type = "list")

  # Tests if it works when CIs are not wanted
  GroupNoCI <- sfmfep(formula = t.formula, method = method, group = "gr", bootstrap = boot, B = B,
                       data = test.data, mu = mu, sigmaCI = NULL,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
  expect_type (object = GroupNoCI, type = "list")

  # Tests unbalanced panels with groups without CI

  test.data <- test.data[-60, ]
  unbalancedGroup <- sfmfep(formula = t.formula, method = method, group ="gr", bootstrap = boot, B = B,
                       data = test.data, mu = mu, sigmaCI = NULL,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
  expect_type (object = unbalancedGroup, type = "list")

  unbalancedBoot <- sfmfep(formula = t.formula, method = method, group ="gr", bootstrap = T, B = 5,
                            data = test.data, mu = mu, sigmaCI = NULL,
                            myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
  expect_type (object = unbalancedBoot, type = "list")

  unbalancedBootNT <- sfmfep(formula = t.formula, method = method, N = 2, Time =c(30,29), bootstrap = T, B = 5,
                           data = test.data, mu = mu, sigmaCI = 0.05,
                           myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
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

context ("SFM.within / SFM.firstDiff")

# Dont really need them anymore, as we already have it tested in the main test function
 test_that ("SFM.within / SFM.firstDiff", {
   sigma_u <- 0.1
   sigma_v <- 0.1
   beta <- c(0.5, 2)
   delta <- c(0.5, 3)

   output <- SFM.within(par = c(sigma_u = sigma_u, sigma_v = sigma_v, beta = beta, delta = delta),
                        xv <- sfm.data[2:3], y <- sfm.data[4], z <- sfm.data[5:6],
                        N <- 2,  Time <- 30, mu = 0, optim = T)  # optim=T for a double (log.l)
   expect_type(object = output, type ="double")

   xv <- matrix(c(rnorm(10,5,2), rnorm(5,20,1), rnorm(10,5,2), rnorm(5,20,1)), ncol=2)
   y <- matrix( c(rnorm(10,5,2), rnorm(5,20,1)))
   z <- matrix(c(rnorm(10,5,2), rnorm(5,20,1), rnorm(10,5,2), rnorm(5,20,1)), ncol=2)
   N <- 2
   Time <- c(10,5)
   data <- data.frame(y = y, xv = xv, z = z)
   output <- SFM.firstDiff(par = c(sigma_u = sigma_u, sigma_v = sigma_v, beta = beta, delta = delta),
                        xv <- xv, y <- y, z <- z,
                        N <- N,  Time <- Time, mu = 0, optim = F)  # optim = F for a list
   expect_type(object = output, type ="list")

   rm(list=ls(all=TRUE))
})
