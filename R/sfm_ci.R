#' < Computes Confidence Intervals for the regression coefficients based on the Hessian Matrix>
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < Describe what is returned when applying this functinon >
#' @examples
#' < an example >


SFM.CI <- function(estimates, hessianMatrix, alpha){

  # Error Handling & exclusion of improper Fisher entries  ---------------------------
  # If the Hessian Matrix is indefinite, we can not calculate Confidence Intervals.
  # Can occure when eigenvalues of the Hessian are != 0 (estimates are saddle points)
  try (fisher_info <- solve(hessianMatrix), silent = T)

  # If diagonal entries of the fisher are negative, these variables are excluded and
  # CIs are calculated for the remaining variables

  try( if (any (diag (fisher_info) < 0) == T){
    indexExcludeVar <- which (diag (fisher_info) < 0)
    fisher_info <- fisher_info[-indexExcludeVar, -indexExcludeVar]
    estimates <- estimates[-indexExcludeVar]
    cat ("Could not compute Confidence Interval for:",
         names(estimates)[indexExcludeVar],"(negative Fisher Information )")
  }, silent = T)

  if (exists ("fisher_info")){  # if dim == 0 no CI is calculated

    standerror <- sqrt( diag (fisher_info))

    # Calculation of CIs  ---------------------------
    # Note that alpha can be a vector of different significance values
    upper <- matrix (c(rep (NA, length (standerror) * length (alpha))), ncol = length (alpha))
    lower <- matrix (c(rep (NA, length (standerror) * length (alpha))), ncol = length (alpha))

    colnames(upper) <- c(1 - alpha[c(1:length (alpha))] / 2)
    colnames(lower) <- c(alpha[c(1:length (alpha))] / 2)

    for (i in 1:length (standerror) && dim (fisher_info)[1] > 0){
      upper[i, ] <- estimates[i] + qnorm (1 - alpha / 2)  * standerror[i]
      lower[i, ] <- estimates[i] + qnorm (alpha / 2) * standerror[i]
    }

    # Create output data.frame for CIs  ---------------------------
    interval <- data.frame (value = estimates, lower = lower, upper = upper, standerror = standerror)
    return (interval)
  } else {
    # User is informed that the Hessian Matrix is not valid.
    # It does not stop() terminate other functions applying this function.
    warning ("Hessian Matrix is singular / indefinite. Could not calculate CIs")
  }
}
