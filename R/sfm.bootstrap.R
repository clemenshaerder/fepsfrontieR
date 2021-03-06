#' @title SFM.bootstrap performs estimation with B Individual Bootstrap Samples
#'
#' @description B Individual Bootstrap Samples are generated from the input and MLE is performed
#'     for each sample. Unlike i.i.d. bootstrapping, individual bootrapping samples
#'     the rows with replacement individually for each panel instead from all panels.
#'     In addition to the mean and the standard error of the estimates,
#'     a confidence interval is returned based on the quantiles of the distribution of estimates.
#' @param y is a n*t x 1 vector (response)
#' @param xv is a n*t x k matrix (explantatory variables)
#' @param z is a n*t x r matrix (inefficency determinants)
#' @param mu is an integer (mean of the truncated normal distribution of the inefficency)
#' @param N is an integer (n - panels)
#' @param Time is an integer (observations per panel)
#' @param method a required string specifying the method ("within" or "firstdiff").
#' @param R is an integer (# of z variables)
#' @param K is an integer (# of xv variables)
#' @param B is an integer (# of Bootstraps)
#' @param myPar is a vecor which has to be entered in the following order:
#'      c(sigma_v, sigma_u, beta = c(), delta = c()).
#'      Required as starting point for the estimation.
#' @param lowerInt is a vector of doubles (lower bound for the estimation)
#' @param alphaCI is a vector of doubles (significance of the Confidence Intervals)
#' @param cumTime ia a vector of the cumulated times of the Time vector.
#'     It serves as an index for computation.
#' @param parallel is an optional boolean variable. If it is set to TRUE, bootstrapping is
#'     performed with parallelization, using all available cores - 1.
#'     Only available for OS Windows.
#' @return A B x (K + R + 2) matrix is returned of the estimates, the mean, standard error
#'     and a confidence interval for each estimate as a data frame.
#' @importFrom parallel detectCores parLapply clusterExport clusterEvalQ stopCluster makeCluster
#' @importFrom dplyr sample_n
#' @importFrom stats quantile


SFM.bootstrap <- function(y, xv, z, mu, N, Time, method, R, K, B,
                          myPar = NULL, lowerInt, alphaCI, cumTime, parallel){

  # create data frame of input variables which helps
  # to conduct the rowise bootstrapping
  data <- data.frame (y = y, xv = xv, z = z)
  rows <- dim (data)[1]
  cols <- dim (data)[2]

  # create the index to sample from the different panels
  index <- findInterval (seq (1:rows), cumTime, left.open = TRUE)

  # draw B individual bootstrap samples. each list entry consists of N lists.
  bootList <- replicate(B, list(), simplify = F)  # create B bootstrap samples

  # for every entry of bootList we sample rowwise for each panel
  bootList <- lapply (bootList, function(x)
                                by (data, simplify = F, INDICES = index,
                                FUN = function(x) sample_n (tbl = x, size = dim (x)[1], replace = T)))

  # transforms the N-lists in the list to a matrices which can be used for bootstrapping
  bootListMat <- lapply(bootList, function(x) do.call (rbind, x))


  if (parallel == F){
    if (method == "within"){
      bootEstimates <- lapply (bootListMat, function(x) nlminb (lower = lowerInt,
                                                                start = myPar,
                                                                Time = Time,
                                                                N = N,
                                                                xv = as.matrix (x[, 2:(2+K-1)]),
                                                                y = as.matrix (x[, 1]),
                                                                z = as.matrix (x[, (2+K):cols]),
                                                                mu = mu,
                                                                optim = T,
                                                                K = K,
                                                                R = R,
                                                                objective = SFM.within,
                                                                # we want only the estimates $par
                                                                cumTime = cumTime)$par)
    } else {
      bootEstimates <- lapply (bootListMat, function(x) nlminb (lower = lowerInt,
                                                                start = myPar,
                                                                Time = Time,
                                                                N = N,
                                                                xv = as.matrix (x[, 2:(2+K-1)]),
                                                                y = as.matrix (x[, 1]),
                                                                z = as.matrix (x[, (2+K):cols]),
                                                                mu = mu,
                                                                optim = T,
                                                                K = K,
                                                                R = R,
                                                                objective = SFM.firstDiff,
                                                                # we want only the estimates $par
                                                                cumTime = cumTime)$par)
    }
  } else {  # parallel = T
      if ((Sys.info()[1] == "Windows") == F){
        stop ("Parallel computing is currently only available for OS Windows.")
      }

      no_of_cores = detectCores()

      # 1 core should have capacity to continue working while bootstrapping
      if(no_of_cores > 1){
        no_of_cores <- no_of_cores - 1
      }

      cl = makeCluster(no_of_cores, type = "PSOCK")

      # We provide each cluster all variables
      optim <- T
      myPar <- myPar
      lowerInt <- lowerInt
      Time <- Time
      N <- N
      bootListMat <- bootListMat
      mu <- mu
      K <- K
      R <- R
      method <- method
      cumTime <- cumTime
      cols <- cols

      # export all variables to the clusters
      clusterExport(cl, c("myPar", "lowerInt", "Time", "N", "cols",
                          "bootListMat", "mu", "optim", "K", "R", "method", "cumTime"),
                    envir = environment())

      # we assign each cluster all required librarys
      clusterEvalQ(cl, {library(fepsfrontieR)})

      if (method == "within"){
        bootEstimates <- parLapply (cl = cl, bootListMat, function(x) nlminb(lower = lowerInt,
                                                                             start = myPar,
                                                                             Time = Time,
                                                                             N = N,
                                                                             xv = as.matrix (x[, 2:(2+K-1)]),
                                                                             y = as.matrix (x[, 1]),
                                                                             z = as.matrix (x[, (2+K):cols]),
                                                                             mu = mu,
                                                                             optim = T,
                                                                             K = K,
                                                                             R = R,
                                                                             objective = SFM.within,
                                                                             # we want only the estimates $par
                                                                             cumTime = cumTime)$par)
      } else {   # use first difference
        bootEstimates <- parLapply (cl = cl, bootListMat, function(x) nlminb(lower = lowerInt,
                                                                             start = myPar,
                                                                             Time = Time,
                                                                             N = N,
                                                                             xv = as.matrix (x[, 2:(2+K-1)]),
                                                                             y = as.matrix (x[, 1]),
                                                                             z = as.matrix (x[, (2+K):cols]),
                                                                             mu = mu,
                                                                             optim = T,
                                                                             K = K,
                                                                             R = R,
                                                                             objective = SFM.firstDiff,
                                                                             # we want only the estimates $par
                                                                             cumTime = cumTime)$par)
      }
      stopCluster (cl)
  }

  # creates a matrix of estimates to calculate colmeans & standard error
  estimatesMat <- do.call (rbind, bootEstimates)

  estimates <- apply (estimatesMat, 2, mean)
  stderror <- apply (estimatesMat, 2, sd)

  # Calculate CIs based on the quantiles of the estimate distribution
  if (!is.null (alphaCI)){
    conf.Interval <- t (apply (estimatesMat, 2,
                               function(x) quantile (x,probs = c(alphaCI/2, 1-alphaCI/2))))
  } else{
    conf.Interval <- NULL
  }

  return(list (estimatesMat = estimatesMat, par = estimates,
               standerror = stderror, conf.Interval = conf.Interval))
}
