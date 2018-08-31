#' @title Inefficiency Index
#' @description Estimator for the inefficiency based on the conditional expectation of "u"
#' on "epsilon". The estimator does not require alpha and thus not suffer from the
#' approximation problem.
#' @param par is a vector of regression coefficients & variance parameters.
#' 1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return returns the mean ineffieciency of each panel


SFM.inindex <- function(h, sigma2star, mu2star, N, Time, method, cumTime){

  #since we only need the square root of sigma two star
  sigmaStar <- sqrt (sigma2star)
  #make sure that operations work
  h <- as.matrix(h)

  # calculates the mean inefficency per panel
  inefficencyIndex <-
    lapply (1:N, function(x)
      mean (h[(cumTime[x] + 1):cumTime[x + 1]] *
              (mu2star[x] + (dnorm (mu2star[x] / sigmaStar[x]) *
                                sigmaStar[x]) / pnorm(mu2star[x] / sigmaStar[x])
              )))

  inefficencyIndex <- as.matrix (unlist (inefficencyIndex))

  return (inefficencyIndex)
}


