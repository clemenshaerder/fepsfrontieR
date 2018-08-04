#' SFM.generate creates a fixed-effects panel stochastic frontier model data set
#' which can be used as examples for other functions of the package fepsfrontieR.
#' @param N is an integer and specifies the amount of panels
#' @param Time is an integer and specifies the amount of observations per panel
#' @param beta is a vector of k estimates.
#' K explenatory varibles will be gerenated.
#' @param delta is a vector of r estimates.
#' R inefficency determinats will be gerenated.
#' @param sigma_u is postivie numeric and is the variation
#' of the stochastic inefficency
#' @param sigma_v is postivie numeric and is the variation
#' of the zero-mean random error.
#' @param mu_u is numeric and is the mean of the truncated normal distribution
#' of the stochastic inefficency.
#' @return A data.fram() including x, y, z & alpha variables
#' @examples
#' < an example >
#' @export

SFM.generate <- function(N, Time, beta, delta, sigma_u, sigma_v, mu_u = 0){

  if (!is.double(N) | !is.double(Time) | !is.vector(beta) | !is.vector(delta) |
      !is.numeric(sigma_u) | !is.numeric(sigma_v)  | !is.numeric(mu_u) |
      sigma_u <= 0 | sigma_v <= 0){
    stop ("Invalid input format of parameters.")
  }

  K <- length(beta)
  R <- length(delta)

  # Generate inefficencys for each panel and repeat it from truncated normal distribution (x>a)
  u_star <- rep(truncnorm::rtruncnorm(N, a = 0, mean = mu_u, sd = sigma_u), each = Time)

  # Generate one alpha intercepts for each panel and repeat it from a uniform distribution
  alpha  <- rep(runif(N, 0, 1), each = Time)

  # Generate data of the model ---------------------------
  v <- c(rnorm(Time*N, 0, sigma_v))
  z <- matrix(rnorm(Time*N*R, 0, 1), nrow = Time*N, ncol = R)  # R inefficency determinants
  x <- matrix(c(rnorm(Time*N*K, rep(alpha, each = K), 1)),
              nrow = Time*N, ncol = K)  # K explenatory variables using alpha as mean

  h <- exp(z %*% delta)
  u <- h * u_star
  epsilon <- v - u

  # Calculation of the response ---------------------------
  y <- alpha + x%*%beta + epsilon

  returnTibble <- dplyr::as_tibble (data.frame (x = x, y = y, z = z, alpha = alpha))
  return (returnTibble)
}
