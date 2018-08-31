#' Stochastic Frontier Model of the GDP of countries
#'
#' @format A data frame with 2296 observations on 7 variables:
#' \describe{
#'   \item{country}{82 countries with observations each}
#'   \item{code}{Each country has its own code}
#'   \item{yr}{The year of the obseration, starting from 1 to 28}
#'   \item{y}{GDP of the countries}
#'   \item{k}{Aggregate physical capital stock}
#'   \item{l}{Labor, the number of individuals in the workforce between
#'       the ages of 15 and 64.}
#'   \item{h}{Human capital adjusted labor, labor force weighted by the mean years of schoolin}
#' }
"panelgdp"
