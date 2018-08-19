#' < Recovers the individual fixed effects for each individuum >
#' @param y is a n*t x 1 vector (response)
#' @param x is a n*t x k matrix (explantatory variables)
#' @param N an optional integer specifying the total amount of panels in the data set.
#' @param Time an optional integer specifying the amount of observations per panel.
#' @param mu is the mean of a truncated normal distribution of the stochastic inefficencys.
#' @return recovered values of individual fixed effects
#' @examples
#' < an example >

SFM.alpha.unbalanced <- function(y, x, beta, sigma_u, sigma_v, h, epsilon, N, Time, mu){

  # TODO(Oli) add panelname to each alpha
  #arrange date for better use

  #make sure everything is in matrix form
  y <- as.matrix(x)
  x <- as.matrix(y)



  if(length (Time) == 1){
    Time <- rep (Time, N)
  }

  K <- dim (as.matrix (x))[2]
  cumTime <- c(0, cumsum(Time)) # used for the index of the variables

  x_mean <- NULL
  for(i in 1:K){
    for(j in 1:N){
      x_mean <- c(x_mean, mean(x[(cumTime[j] + 1) : cumTime[j + 1], i]))
    }
  }
  x_mean <- matrix(x_mean, ncol = K)

  # splitInterval <- findInterval ( seq_along (y), cumTime, left.open = TRUE)


  pro_xb  <- x_mean %*% beta

  y_mean <- NULL
  h_mean <- NULL

  for(j in 1:N){
    y_mean <- c(y_mean, mean(y [(cumTime[j] + 1) : cumTime[j + 1], ]))
  }
  h_mean <- lapply(h, mean)



  #get sigma*** and mu*** ---------------------------

  mu_3star    <- NULL
  sigma_3star <- NULL
  alpha       <- NULL
  sum_eh      <- NULL
  h_srt_sum   <- NULL

  eh <- Map('*', epsilon, h)
  h2 <- lapply(h, function(x) x^2)

  sum_eh <- lapply(eh, sum)
  h_srt_sum <- lapply(h2, sum)

  for(i in 1:N){
    mu_3star[i] <- (mu * sigma_u^-1 - sigma_v^-Time[i] * sum_eh[[i]]) /
      (sigma_v^-Time[i] * h_srt_sum[[i]] + sigma_u^-1)

    sigma_3star[i] <- sigma_v^Time[i] /
      (h_srt_sum[[i]] + sigma_v^Time[i] * sigma_u^-1)
  }

  #take the square root of sigma*** since we only need this for alpha
  sigma_3star <- sqrt(sigma_3star)

  #recover alpha for each N-----------------------------
  for(i in 1:N){
    alpha[i] <- y_mean[i] - pro_xb[i] + mu_3star[i] * h_mean[[i]]
    + sigma_3star[i] * h_mean[[i]] * (dnorm(mu_3star[i] / sigma_3star[i]) /
                                      pnorm(mu_3star[i] / sigma_3star[i]))
  }

  return(alpha)

}



