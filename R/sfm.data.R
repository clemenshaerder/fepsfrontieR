#' Randomly generated fixed-effects panel stochastic frontier model data set
#'
#' The data set was created using the in-package function SFM.generate with
#' parameters:
#' N = 5, Time = 10,
#' beta <- c(0.5); K <- length(beta);
#' delta <- c(0.5); R <- length(delta)
#' sigma_u <- 0.2; sigma_v <- 0.1
#' set.seed(665)
#'
#'
#' @format A data frame with 60 rows and k+r+1 variables:
#' \describe{
#'   \item{producer}{ID of 10 different producers}
#'   \item{x1}{First explenatory variable}
#'   \item{x2}{Second explenatory variable}
#'   \item{y}{The response variable}
#'   \item{z1}{First inefficency determinant}
#'   \item{z2}{Second inefficency determinant}
#' }
"sfm.data"
