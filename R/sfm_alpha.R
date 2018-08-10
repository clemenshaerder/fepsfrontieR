#' < Recovers the alpha (inteysfdysyxyssfdrcept) of each panel >
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < Describe what is returned when applying this functinon >
#' @examples
#' < an example >


SFM.alpha <- function(y, x, beta, sigma_u, sigma_v, h, epsilon, N, Time, mu = 0){

  # TODO(Oli) comment the file according googles r comment standard (use -------- to easily readable code bloks)
  # TODO(Oli) add panelname to each alpha
  ####-------------------####

  epsilon <- matrix(epsilon, ncol = N , nrow = Time)

  h <- matrix(h, ncol = N, nrow = Time)
  h_mean <- apply(h, 2, mean)

  y <- matrix(y, ncol = N, nrow = Time)
  y_mean <- apply(y, 2, mean)

  x.mat <- matrix(x, nrow = Time)
  x_MEAN <- apply(x.mat, 2, mean)
  x_mean <- matrix(x_MEAN, nrow = N)
  pro_xb <- x_mean %*% beta

  ####--------------------####

  mu_3star <- NULL; sigma_3star <- NULL; alpha <- NULL

  sum_eh <- apply(epsilon*h, 2, sum)
  h_srt_sum <- apply(h^2, 2, sum)

  for(i in 1:N){
    mu_3star[i] <- (mu * sigma_u^-1 - sigma_v^-Time * sum_eh[i]) /
      (sigma_v^-Time * h_srt_sum[i] + sigma_u^-1)

    sigma_3star[i] <- sigma_v^Time /
      (h_srt_sum[i] + sigma_v^Time * sigma_u^-1)
  }

  sigma_3star <- sqrt(sigma_3star)

  for(i in 1:N){
    alpha[i] <- y_mean[i] - pro_xb[i] + mu_3star[i] * h_mean[i]
    + sigma_3star[i] * h_mean[i] * (dnorm(mu_3star[i] / sigma_3star[i]) /
                                      pnorm(mu_3star[i] / sigma_3star[i]))
  }

  return(alpha)

}
