#' @title Inefficiency Index of the Stochastic Frontier Model
#' @description Estimation of the inefficiency based on the conditional expectation of "u"
#'     on "epsilon". The estimator does not require alpha and thus not suffer from the
#'     approximation problem.
#' @param h are the values of the positive non-stochastic inefficency determinants.
#' @param sigma2star is a vector of variance estimations based on a chosen transformation.
#' @param mu2star is a vector of mean estimations based on a chosen transformation.
#' @param N an optional integer specifying the total amount of panels in the data set.
#' @param cumTime ia a vector of the cumulated times of the Time vector.
#'     It serves as an index for computation.
#' @return returns the mean ineffieciency of each panel as a Nx1 vector


SFM.inindex <- function(h, sigma2star, mu2star, N, cumTime){

  #since we only need the square root of sigma two star
  sigmaStar <- sqrt (sigma2star)
  #make sure that operations work
  h <- as.matrix(h)

  # calculates the mean inefficency per panel
  inefficencyIndex <- lapply (1:N, function(x)
                                   mean (h[(cumTime[x] + 1):cumTime[x + 1]] *
                                   (mu2star[x] + (dnorm (mu2star[x] / sigmaStar[x]) *
                                   sigmaStar[x]) / pnorm(mu2star[x] / sigmaStar[x]))))

  inefficencyIndex <- as.matrix (unlist (inefficencyIndex))

  return (inefficencyIndex)
}


