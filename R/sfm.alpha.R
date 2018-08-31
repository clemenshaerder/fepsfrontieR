#' @title Alpha
#' @description Recovers the individual fixed effects by using the first-order condition for
#' alpha from the untransformed log-likelihood function of the model assuming
#' all other parameters are known. Altough, the variance of alpha is high for small
#' sample sized date.
#' @param y is a n*t x 1 vector (response)
#' @param x is a n*t x k matrix (explantatory variables)
#' @param N an optional integer specifying the total amount of panels in the data set.
#' @param Time an optional integer specifying the amount of observations per panel.
#' @param mu is the mean of a truncated normal distribution of the stochastic inefficencys.
#' @return recovered values of individual fixed effects


SFM.alpha <- function(y, x, beta, sigma_u, sigma_v, h, epsilon, N, Time, mu, cumTime){

  # dimensions of x
  K <- dim (as.matrix (x))[2]

  # Calculate the mean over all overservations within each panel of x
  x_mean <- NULL
  for (i in 1:K) {
    for (j in 1:N) {
      x_mean <- c(x_mean, mean (x[(cumTime[j] + 1):cumTime[j + 1], i]))
    }
  }
  x_mean <- matrix (x_mean, ncol = K)

  # splitInterval <- findInterval ( seq_along (y), cumTime, left.open = TRUE)

  # Calculate the mean over all observations within each panel of y
  y_mean <- c()
  y_mean <- lapply (1:N, function(x)
    y_mean <-
      c(y_mean, mean (y [(cumTime[x] + 1):cumTime[x + 1],])))

  y_mean <- as.matrix (unlist (y_mean))


  # calculate sigma*** and mu*** for each panel ---------------------------

  # mean of h of each panel
  h_mean    <- lapply (h, mean)
  # sum of each panels of h where each element in the panel is times 2
  h2_sum    <- lapply (h, function(x) sum(x ^ 2))
  # multiplies the each element of epsilon and h of the same panel
  # and build the sum of it
  eh        <- Map ('*', epsilon, h)
  sum_eh    <- lapply (eh, sum)

  #get mu_3star
  mu_3star <-
    lapply (1:N, function(x)
      (mu * sigma_u ^ -1 - sigma_v ^ -Time[x] * sum_eh[[x]]) /
        (sigma_v ^ -Time[x] * h2_sum[[x]] + sigma_u ^ -1))
  mu_3star <- as.matrix(unname (unlist (mu_3star)))

  #get sigma_3star
  sigma_3star <-
    lapply (1:N, function(x)
      sigma_v ^ Time[x] / (h2_sum[[x]] + sigma_v ^
                             Time[x] * sigma_u ^ -1))
  sigma_3star <- as.matrix(unname (unlist (sigma_3star)))

  # take the square root of sigma*** since we only need this for alpha
  sqrt_sigma_3star <- sqrt (sigma_3star)


  # recover alpha for each panel N -----------------------------

  pro_xb  <- x_mean %*% beta

  #main operation to recover alpha from estimates
  alpha <-
    lapply (1:N, function(x)
      y_mean[x] - pro_xb[x] + mu_3star[x] * h_mean[[x]] +
        sqrt_sigma_3star[x] * h_mean[[x]] *
        if (sqrt_sigma_3star[x] == 0 ||
            pnorm (mu_3star[x] / sqrt_sigma_3star[x]) == 0) {
          1
        } else {
          (dnorm (mu_3star[x] / sqrt_sigma_3star[x]) /
             pnorm (mu_3star[x] / sqrt_sigma_3star[x]))

        })

  alpha <- as.matrix (unname(unlist(alpha)))

  return (alpha)
}

