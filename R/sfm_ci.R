#' < Computes Confidence Intervals for the regression coefficients based on the Hessian Matrix>
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < Describe what is returned when applying this functinon >
#' @examples
#' < an example >


SFM.CI <- function(estimates, hessianMatrix, alpha){

  # Error Handling: if the Hessian Matrix is indefinite, we can not calculate Confidence Intervals
  # Can occure when eigenvalues of the Hessian are != 0 (estimates are saddle points)
  try (fisher_info <- solve(hessianMatrix), silent = T)
  try (prop_sigma <- sqrt(diag(fisher_info)), silent=T)  # checks if variance is improper (-)

  if (exists("prop_sigma")){

    # Calculation of CIs  ---------------------------
    # Note that alpha can be a vector of different significance values
    upper <- matrix(c(rep(NA, length(estimates)* length(alpha))), ncol = length(alpha))
    lower <- matrix(c(rep(NA, length(estimates)* length(alpha))), ncol = length(alpha))

    colnames(upper) <- c(1-alpha[c(1:length(alpha))]/2)
    colnames(lower) <- c(alpha[c(1:length(alpha))]/2)

    for (i in 1:length(estimates)){
      upper[i, ] <- estimates[i]+qnorm(1-alpha/2)*prop_sigma[i]
      lower[i, ] <- estimates[i]+qnorm(alpha/2)*prop_sigma[i]
    }

    # Create output data.frame for CIs  ---------------------------
    interval <- data.frame(value = estimates, lower = lower, upper = upper)
    return (interval)
  } else {
    # Error handling of SFM.CI  ---------------------------
    # User is informed that the Hessian Matrix is not valid.
    # It does not stop() terminate other functions applying this function.
    warning ("Hessian Matrix is singular / indefinite. Could not calculate CIs")
  }

}
