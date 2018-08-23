#' SFM.bootstrap performs estimation with B Individual Bootstrap Samples
#'
#' B Individual Bootstrap Samples are genereated from the input and MLE is performed
#' for each sample. Unlike i.i.d. bootstrapping, individual bootrapping samples
#' the rows with replacement individually for each panel instead from all samples.
#' In addition to the mean and the standard error of the estimates,
#' a confidence interval is returned based on the quantiles of the distribution of estimates.
#' @param B is an integer (# of Bootstraps)
#' @param xv is a n*t x k matrix (explantatory variables)
#' @param z is a n*t x r matrix (inefficency determinants)
#' @param y is a n*t x 1 vector (response)
#' @param mu is an integer (mean of the truncated normal distribution of the inefficency)
#' @param R is an integer (# of z variables)
#' @param K is an integer (# of xv variables)
#' @param N is an integer (n - panels)
#' @param Time is an integer (observations per panel)
#' @param method a required string specifying the method ("within" or "firstdiff").
#' @param lowerInt is a vector of doubles (lower bound for the estimation)
#' @param sigmaCI is a vector of doubles (significance of the Confidence Intervals)
#' @param myPar is a vecor which has to be entered in the following order:
#'      c(sigma_v, sigma_u, beta = c(), delta = c()). Required as starting point for the estimation.
#' @return A B x (K + R + 2) matrix is returned of the estimates, the mean, standard error
#'      and a confidence interval for each estimate as a data frame.

# TODO: estimates are off. try to recreate bootstrapping

SFM.bootstrap <- function(y, xv, z, mu, N, Time, method, R, K, B, myPar = NULL, lowerInt, sigmaCI, cumTime){

  # create data frame of input variables which helps
  # to conduct the rowise bootstrapping
  data <- data.frame (y = y, xv = xv, z = z)
  rows <- dim (data)[1]
  cols <- dim (data)[2]

  # create the index to sample from the different panels
  index <- findInterval (seq (1:rows), cumTime, left.open = TRUE)

  # draw R individual bootstrap samples. each list entry consists of N lists.
  bootList <- replicate(B, list(), simplify = F)  # create R bootstrap samples
  # for every entry of bootList we sample rowwise for each panel
  bootList <- lapply (bootList, function(x)
                                by (data, simplify = F, INDICES = index,
                                FUN = function(x) sample_n (tbl = x, size = dim (x)[1], replace = T)))

  # transforms the N-lists in the list to a matrices which can be used for bootstrapping
  bootListMat <- lapply(bootList, function(x) do.call (rbind, x))

  # TODO(): OPTIMIZATION by parallization??
  # parallel::parLapply(lapply (bootList, function(x) by (data, INDICES = index, FUN = function(x) sample_n (tbl = x, size = dim(x)[1], replace = T))))

  # calculate estimates for each entry of the list depending on the method
  if (method == "within"){
    bootEstimates <- lapply (bootListMat, function(x) nlminb(lower = lowerInt,
                                                      start = myPar,  # TBD by Rouven
                                                      Time = Time,
                                                      N = N,
                                                      xv = as.matrix (x[, 2:(2+K-1)]),
                                                      y = as.matrix (x[, 1]),
                                                      z = as.matrix (x[, (2+K):cols]),
                                                      mu = mu,
                                                      optim = T,
                                                      K = K, R = R,
                                                      objective = SFM.within,
                                                      cumTime = cumTime
                                                      )$par)  # we want only the estimates
  } else {
    bootEstimates <- lapply (bootListMat, function(x) nlminb(lower = lowerInt,
                                                             start = myPar,  # TBD by Rouven
                                                             Time = Time,
                                                             N = N,
                                                             xv = as.matrix (x[, 2:(2+K-1)]),
                                                             y = as.matrix (x[, 1]),
                                                             z = as.matrix (x[, (2+K):cols]),
                                                             mu = mu,
                                                             optim = T,
                                                             K = K, R = R,
                                                             objective = SFM.firstDiff,
                                                             cumTime = cumTime
                                                             )$par)  # we want only the estimates
  }

  # creates a matrix of estimates to calculate colmeans & standard error
  estimatesMat <- do.call (rbind, bootEstimates)

  estimates <- apply (estimatesMat, 2, mean)
  stderror <- apply (estimatesMat, 2, sd)

  # Calculate CIs based on the quantiles of the estimate distribution
  if (!is.null(sigmaCI)){
    conf.Interval <- t (apply (estimatesMat, 2, function(x) quantile(x, probs = c(sigmaCI/2, 1-sigmaCI/2))))
  } else{
    conf.Interval <- "NULL"
  }

  # TODO() we could include a histogram of the estimates and QQ-Plot.
  # would be nice but not a must.
  return(list (estimatesMat = estimatesMat, par = estimates,
               standerror = stderror, conf.Interval = conf.Interval))
}

