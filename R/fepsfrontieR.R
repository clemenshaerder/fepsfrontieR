#' Estimating Fixed-Effects Panel Stochastic Frontier Models
#'
#' A package for estimation of fixed-effect panel stochastic
#' frontier models by within or first-difference model transformation
#' as proposed by Wang & Ho (2010).
#' Maximum Likelihood Estimation is performed using nlminb( ).
#' Standard Errors are provided using the Hessian Matrix or
#' by individual Bootstrapping.
#' The fepsfrontieR package provides the main function sfmfep
#' and a supplementary data generation function SFM.generate
#'
#' @section sfmfep:
#' The sfmfep functions is the main function of the package fepsfrontieR.
#' sfmfep can be either used to estimate model parameters based on a data input,
#' or to fit a prespecified model.
#' Estimation is performed by applying a within or first-difference transformation
#' of the data. One can optionally chose Bootstrapping to calculate the estimates and
#' standard errors of the estimates. Otherwise the Standard Errors are obtained
#' from the Hessian matrix.
#' The estimation of the parameters is numerically done using the optimizer nlminb ( ).
#' Starting points for the optimizer can be provided. Otherwise, OLS estimators
#' are used as starting points.
#' Finally, the inefficency index & fixed-effects are calculated.
#' A S3 object is returned.
#'
#' @section SFM.generate:
#' The SFM.generate function generates a balanced fixed-effects stochastic
#' frontier model with exponential inefficencys. The amount of specified
#' betas & deltas create associated explenatory varibles & inefficeny determinants.
#'
#' @docType package
#' @name fepsfrontieR
NULL
