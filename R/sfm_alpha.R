#' < Recovers the individual fixed effects for each individuum >
#' @param y is a n*t x 1 vector (response)
#' @param x is a n*t x k matrix (explantatory variables)
#' @param N an optional integer specifying the total amount of panels in the data set.
#' @param Time an optional integer specifying the amount of observations per panel.
#' @param mu is the mean of a truncated normal distribution of the stochastic inefficencys.
#' @return recovered values of individual fixed effects
#' @examples
#' < an example >
#' @export

SFM.alpha <- function(y, x, beta, sigma_u, sigma_v, h, epsilon, N, Time, mu, cumTime){

  K <- dim (as.matrix (x))[2]

  # Calculate the mean of each panel for each explenatory variable
  x_mean <- NULL
  for(i in 1:K){
    for(j in 1:N){
      x_mean <- c(x_mean, mean(x[(cumTime[j] + 1) : cumTime[j + 1], i]))
    }
  }
  x_mean <- matrix(x_mean, ncol = K)

  # splitInterval <- findInterval ( seq_along (y), cumTime, left.open = TRUE)
  y_mean <- c()
  y_mean <- lapply (1:N, function(x)
                         y_mean <- c(y_mean, mean(y [(cumTime[x] + 1) : cumTime[x + 1], ])))
  y_mean <- as.matrix (unlist (y_mean))


  # calculate sigma*** and mu*** for each panel ---------------------------

  h_mean    <- lapply (h, mean)  # mean of h of each panel
  eh        <- Map('*', epsilon, h)  # multiplies the lists epsilon & h
  h2_sum    <- lapply(h, function(x) sum(x^2))
  sum_eh    <- lapply (eh, sum)

  mu_3star <- lapply (1:N, function(x) (mu * sigma_u^-1 - sigma_v^-Time[x] * sum_eh[[x]]) /
                        (sigma_v^-Time[x] * h2_sum[[x]] + sigma_u^-1))
  mu_3star <- as.matrix (unname (unlist (mu_3star)))

  sigma_3star <- lapply (1:N, function(x) sigma_v^Time[x] / (h2_sum[[x]] +
                                                               sigma_v^Time[x] * sigma_u^-1))
  sigma_3star <- as.matrix (unname (unlist (sigma_3star)))

  # take the square root of sigma*** since we only need this for alpha
  sqrt_sigma_3star <- sqrt(sigma_3star)


  # recover alpha for each panel N -----------------------------

  pro_xb  <- x_mean %*% beta

  alpha <- lapply (1:N, function(x) y_mean[x] - pro_xb[x] + mu_3star[x] * h_mean[[x]] +
                     sqrt_sigma_3star[x] * h_mean[[x]] *
                     # if values are 0 we set it to 1  to get a valid result
                     if (sqrt_sigma_3star[x] == 0 ||
                         pnorm (mu_3star[x] / sqrt_sigma_3star[x] ) == 0){
                            1
                       } else {
                       (dnorm (mu_3star[x] / sqrt_sigma_3star[x]) /
                          pnorm (mu_3star[x] / sqrt_sigma_3star[x]))
                       }
                   )

  alpha <- as.matrix (unname (unlist (alpha)))

  return(alpha)

}

