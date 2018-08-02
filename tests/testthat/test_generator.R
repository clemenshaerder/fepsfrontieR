context ("Correct output of SFM.generate")

library (MASS)
library (dplyr)

test_that ("sfmfep works :-)", {

  N <- 2
  Time <- 30
  beta <- c(0.5,3)
  delta <- c(0.5,2)
  sigma_u <- 2; sigma_v <- 1
  form.test <- formula(y  ~ x1 + x2 + (z1 + z2))
  test.data <- sfm.data

  testSfmfep <- sfmfep(formula = form.test, N=2,Time=30, data = test.data, mu = 0, optimPar = NULL)
  expect_type (object = testSfmfep, type = "list")
})


context ("Correct output of SFM.generate")
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
                           mu_u = mu, sigma_u = sigma_u, sigma_v = sigma_v)
  expect_type (object = testSFM1, type = "list")

  N <- "c(3,2)"
  expect_error( SFM.generate(N = N, Time = Time, beta = beta, delta = delta,
                             mu_u = mu, sigma_u = sigma_u, sigma_v = sigma_v))

  rm(list=ls(all=TRUE))
})




