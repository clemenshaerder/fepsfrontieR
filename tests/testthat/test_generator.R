context ("sfmfep")

library (MASS)
library (dplyr)

test_that ("sfmfep works", {

  N <- 2
  Time <- 30
  beta <- c(0.5,3)
  delta <- c(0.5,2)
  sigma_u <- 2; sigma_v <- 1
  form.test <- formula(y  ~ x1 + x2 + (z1 + z2))
  test.data <- sfm.data  # package data

  testSfmfep <- sfmfep(formula = form.test, method= "within", N=2,Time=30, data = test.data, mu = 0, myPar = NULL)
  expect_type (object = testSfmfep, type = "list")

  # Tests if optim " N & T" works
  testSfmfep <- sfmfep(formula = form.test, N=2,Time=30, method = "firstdiff", data = test.data, mu = 0, myPar = NULL)
  expect_type (object = testSfmfep, type = "list")

  # Tests if optim " N & T as unbalanced vector" works
  expect_error ( sfmfep(formula = form.test, method = "firstdiff", N=2,Time=c(30,29), data = test.data, mu = 0, myPar = NULL) )

  # Tests if option "group" TODO:(currently throws a NaN but result is correct)
  testSfmfep <- sfmfep(formula = form.test, method = "firstdiff", group ="gr", data = test.data, mu = 0, myPar = NULL)
  expect_type (object = testSfmfep, type = "list")

  # Tests if defined starting points "myPar" works
  testSfmfep <- sfmfep(formula = form.test, method = "firstdiff", group ="gr",
                       data = test.data, mu = 0,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
  expect_type (object = testSfmfep, type = "list")

  # Tests if it works when estimates provided (estimate = F)
  testSfmfep <- sfmfep(formula = form.test, method = "firstdiff", group ="gr",
                       data = test.data, mu = 0, estimate = F,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
  expect_type(object = testSfmfep, type = "list")

  # Tests if it works when CIs are not wanted
  testSfmfep <- sfmfep(formula = form.test, method = "firstdiff", group ="gr",
                       data = test.data, mu = 0, sigmaCI = NULL,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
  expect_type (object = testSfmfep, type = "list")

  # Tests unbalanced panels
  test.data <- test.data[-60, ]
  testSfmfep <- sfmfep(formula = form.test, method = "firstdiff", group ="gr",
                       data = test.data, mu = 0, sigmaCI = NULL,
                       myPar = c(sigma_u = 1, sigma_v=2, beta = c(1,2), delta = c(1, 2)))
  expect_type (object = testSfmfep, type = "list")
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

context ("SFM.within")

 test_that ("SFM.within / SFM.firstDiff", {
   sigma_u <- 0.1
   sigma_v <- 0.1
   beta <- c(0.5, 2)
   delta <- c(0.5, 3)

   output <- SFM.firstDiff(par = c(sigma_u = sigma_u, sigma_v = sigma_v, beta = beta, delta = delta),
                        xv <- sfm.data[2:3], y <- sfm.data[4], z <- sfm.data[5:6],
                        N <- 2,  Time <- 30, mu = 0, optim = T)  # optim=T for a double (log.l)
   expect_type(object = output, type ="double")

   xv <- matrix(c(rnorm(10,5,2), rnorm(5,20,1), rnorm(10,5,2), rnorm(5,20,1)), ncol=2)
   y <- matrix( c(rnorm(10,5,2), rnorm(5,20,1)))
   z <- matrix(c(rnorm(10,5,2), rnorm(5,20,1), rnorm(10,5,2), rnorm(5,20,1)), ncol=2)
   N <- 2
   Time <- c(10,5)
   output <- SFM.firstDiff(par = c(sigma_u = sigma_u, sigma_v = sigma_v, beta = beta, delta = delta),
                        xv <- xv, y <- y, z <- z,
                        N <- N,  Time <- Time, mu = 0, optim = F)  # optim = F for a list
   expect_type(object = output, type ="list")


})
