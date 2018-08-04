#' < Recovers the alpha (inteysfdysyxyssfdrcept) of each panel >
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < Describe what is returned when applying this functinon >
#' @examples
#' < an example >

SFM.alpha <- function(y, x, beta, sigma_u, sigma_v, h, epsilon, N, Time, mu = 0){

  # TODO(Oli) comment the file according googles r comment standard (use -------- to easily readable code bloks)
  # TODO(Oli) add panelname to each alpha

  sigma_u <- sqrt(sigma_u) # since we need sigma times - 2 and i don`t know how to get it in another way`
  sigma_v <- sqrt(2)  # TODO(Oli) warum 2? sigma_v?

  # converge h into a matrix and get the mean for each panel
  h <- matrix(h, ncol = N, nrow = Time)
  h_mean <- apply(h, 2, mean)

  #and also epsilon
  epsilon <- matrix(epsilon, ncol = N , nrow = Time)

  # get the mean of y´s panels
  y <- matrix(y, ncol = N, nrow = Time)
  y_mean <- apply(y, 2, mean)

  # get the mean of xÂ´s panels
  x.mat <- matrix(x, nrow = Time)
  x_mean <- apply(x.mat, 2, mean)
  x_mean <- matrix(x_mean, nrow = N)

  # create a NULL matrix which get filled in the next step
  a <- matrix(rep(NA, N*Time), nrow = Time, ncol = N)

  for (i in 1:N){            # creates a matrix with the product of epsilon and h for each Timepoint for each panel
    for (j in 1:Time){
       a[j, i] <- epsilon[j, i] * h[j, i]
    }
  }

  b <- apply(a, 2, sum)     # creates a vector with the sum of each row of a

  nominator   <- NULL # TODO(Oliver): maybe get rid of nominator & denominator
  denominator <- NULL
  mu3star     <- NULL
  sigma3star  <- NULL #empty vectors which get filled in the loop

  h_srt_sum <- apply(h^2, 2 , sum) # creates the squared sum of each row for all n

  for(i in 1:N){ # TODO(Oliver): maybe get rid of nominator & denominator -> see within model
    nominator[i] <- mu * sigma_u^(-2) - sigma_v^(-2*Time) * b[i] # a vector with every product of our nominator
    denominator[i] <- sigma_v^(-2*Time) * h_srt_sum[i] + sigma_u^(-2)
    mu3star[i] <- nominator[i] / denominator[i]
    sigma3star[i] <- sigma_v^(2*Time) / ((h_srt_sum[i] + sigma_v^(2*Time) * sigma_u^(-2)))
  }

  ############################################################
  alpha <- NULL
  x_b <- x_mean %*% beta
  summary(alpha)
  for(i in 1:N){
    alpha[i] <- y_mean[i] - x_b[i] + mu3star[i] * h_mean[i] + sqrt(sigma3star[i]) *
                h_mean[i] * ((dnorm(mu3star[i]/ sqrt(sigma3star[i]))) /
                (pnorm(mu3star[i]/ sqrt(sigma3star[i]))))
  }

  return(alpha)
}
