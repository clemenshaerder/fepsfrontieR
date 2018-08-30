#' Randomly generated fixed-effects panel stochastic frontier model data set
#'
#' The data set was created using the in-package function SFM.generate with
#' parameters:
#' N = 2, Time = 30,
#' beta <- c(0.5,3); K <- length(beta);
#' delta <- c(0.5,2); R <- length(delta)
#' sigma_u <- 2; sigma_v <- 1
#' set.seed(665)
#'
#'
#' @format A data frame with 60 rows and k+r+1 variables:
#' \describe{
#'   \item{gr}{2 groups of producers (unsorted)}
#'   \item{x1}{First explenatory variable}
#'   \item{x2}{Second explenatory variable}
#'   \item{y}{The response variable}
#'   \item{z1}{First inefficency determinant}
#'   \item{z2}{Second inefficency determinant}
#' }
"sfm.data"
