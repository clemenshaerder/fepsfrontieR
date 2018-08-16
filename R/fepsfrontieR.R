#' Estimating Fixed-Effects Panel Stochastic Frontier Models
#'
#' A package for estimation of fixed-effect panel stochastic
#' frontier models by within-model transformation as proposed by Wang & Ho (2010).
#' Maximum Likelihood Estimation is performed using nlminb( ).
#' Confidence Intervals are provided using the Hessian Matrix or
#' numerically expensive Bootstrapping to calculate Standard Error.
#' In addition, a Gibbs Sampler can be chosen as an alternative to MLE.
#'
#' The fepsfrontieR package provides one important function
#' sfmfep and a supplementary data generation function SFM.generate
#'
#' @section sfmfep:
#' The sfmfep functions is the main function of the package fepsfrontieR.
#' It provides two options.
#' If "estimate = T" is chosen, the model is
#' estimated with MLE or a Gibbs Sampler. Starting points can be provided
#' for the numerical optimizer nlminb ( ). Otherwise OLS estimators are
#' applied as starting points.
#'
#' If "estimate = F", a stochastic frontier model is calculated,
#' based on provided point estimates.
#'
#' A S3 object is returned in both cases.
#'
#' @usage
#' sfmfep(formula, data, group = NULL, N = NULL,
#'        Time = NULL, mu = 0,  sigmaCI = 0.05, estimate = T,
#'        myPar = c(sigma_u, sigma_v, beta = c( ), delta = c( )))
#'
#' SFM.generate(N, Time, beta, delta, sigma_u, sigma_v, mu_u = 0)
#'
#' @section SFM.generate:
#' The SFM.generate function generates a fixed-effects stochastic
#' frontier model with exponential inefficencys. The amount of specified
#' betas & deltas create associated explenatory varibles & inefficeny determinants.
#'
#' @docType package
#' @name fepsfrontieR
