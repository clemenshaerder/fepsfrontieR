#' < calculate the inefficiency of each individuum >
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < returns the inefficiency for each individuum >
#' @examples
#' < an example >

SFM.inindex.unbalanced <- function(h, sigma2star, mu2star, N, Time){

  # TODO(Oli) add panelname to each ineff.index
  sigma <- sqrt(sigma2star)
  cumTime <- c(0, cumsum(Time))
  h <- as.matrix(h)
  in_index <- c()

  for (i in 1:N){
    in_index <- c(in_index, h[(cumTime[i] + 1) : cumTime[i + 1]] *
                 (mu2star[i] + (dnorm(mu2star[i] / sigma[i]) *
                  sigma[i]) / pnorm(mu2star[i]/ sigma[i])))
  }

  return (in_index)

}
library("gamlss")


