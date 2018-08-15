#' < Recovers the individual fixed effects for each individuum >
#' @param y is a n*t x 1 vector (response)
#' @param x is a n*t x k matrix (explantatory variables)
#' @param N an optional integer specifying the total amount of panels in the data set.
#' @param Time an optional integer specifying the amount of observations per panel.
#' @param mu is the mean of a truncated normal distribution of the stochastic inefficencys.
#' @return recovered values of individual fixed effects
#' @examples
#' < an example >


SFM.alpha <- function(y, x, beta, sigma_u, sigma_v, h, epsilon, N, Time, mu){

  # TODO(Oli) add panelname to each alpha
  #arrange date for better use

  epsilon <- matrix (epsilon, ncol = N , nrow = Time)

  h       <- matrix (h, ncol = N, nrow = Time)
  h_mean  <- apply (h, 2, mean, na.rm = T)

  y       <- matrix (y, ncol = N, nrow = Time)
  y_mean  <- apply (y, 2, mean, na.rm = T)

  x.mat   <- matrix (x, nrow = Time)
  x_MEAN  <- apply (x.mat, 2, mean, na.rm = T)
  x_mean  <- matrix (x_MEAN, nrow = N)
  pro_xb  <- x_mean %*% beta


  #get sigma*** and mu*** ---------------------------

  mu_3star    <- NULL
  sigma_3star <- NULL
  alpha       <- NULL

  sum_eh    <- apply (epsilon*h, 2, sum, na.rm = T)
  h_srt_sum <- apply (h^2, 2, sum, na.rm = T)

  for(i in 1:N){
    mu_3star[i] <- (mu * sigma_u^-1 - sigma_v^-Time * sum_eh[i]) /
      (sigma_v^-Time * h_srt_sum[i] + sigma_u^-1)

    sigma_3star[i] <- sigma_v^Time /
      (h_srt_sum[i] + sigma_v^Time * sigma_u^-1)
  }

  #take the square root of sigma*** since we only need this for alpha
  sigma_3star <- sqrt(sigma_3star)

  #recover alpha for each N-----------------------------
  for(i in 1:N){
    alpha[i] <- y_mean[i] - pro_xb[i] + mu_3star[i] * h_mean[i]
    + sigma_3star[i] * h_mean[i] * (dnorm(mu_3star[i] / sigma_3star[i]) /
    pnorm(mu_3star[i] / sigma_3star[i]))
  }

  return(alpha)

}


