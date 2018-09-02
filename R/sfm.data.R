#' Randomly generated fixed-effects panel stochastic frontier model data set
#'
#' The data set was created using the in-package function SFM.generate with
#' parameters:
#' N = 40, Time = 5,
#' beta <- c(0.5); K <- length(beta);
#' delta <- c(0.5); R <- length(delta)
#' sigma_u <- 0.2; sigma_v <- 0.1
#' set.seed(665)
#'
#'
#' @format A data frame with 60 rows and k+r+1 variables:
#' \describe{
#'   \item{producer}{ID of 20 different producers}
#'   \item{x}{An explenatory variable}
#'   \item{y}{The response variable}
#'   \item{z}{An inefficency determinant}
#'   \item{alpha}{The fixed-effects per producer}
#' }
"sfm.data"
