#' < What it does >
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < Describe what is returned when applying this functinon >
#' @examples
#' < an example >

SFM.inindex <- function(h, sigma2star, mu2star, N, Time){

  # TODO(Oli) comment the file according googles r comment standard (use -------- to easily readable code bloks)
  # TODO(Oli) add panelname to each ineff.index

  h <- matrix(h, ncol = N, nrow = Time)
  in_index <- matrix(rep(NA, N*Time), ncol = N, nrow = Time)

  for (i in 1:N){
    in_index[, i] <- h[, i] * (mu2star[i] + (dnorm(mu2star[i] / sqrt(sigma2star[i])) *
                    sqrt(sigma2star[i])) / pnorm(mu2star[i]/sqrt(sigma2star[i])))
  }

  return (as.vector(in_index))

}
