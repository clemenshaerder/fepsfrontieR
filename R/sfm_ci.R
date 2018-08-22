#' < Computes Confidence Intervals for the regression coefficients based on the Hessian Matrix>
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < Describe what is returned when applying this functinon >
#' @examples
#' < an example >
#' @export


SFM.CI <- function(estimates, hessianMatrix, alpha, N, Time, df){

  # Exclusion of improper Fisher entries  ---------------------------

  # If the Hessian Matrix is indefinite, we can not calculate Confidence Intervals.
  # Can occure when eigenvalues of the Hessian are != 0 (estimates are saddle points)
  try (fisher_info <- solve(hessianMatrix), silent = T)

  # If diagonal entries of the fisher are negative, these variables are excluded and
  # CIs are calculated for the remaining variables

  if (any (diag (fisher_info) < 0) == T){

    indexExcludeVar <- which (diag (fisher_info) < 0)
    fisher_info <- fisher_info[-indexExcludeVar, -indexExcludeVar]

    cat ("Could not compute Confidence Interval for:",
         names(estimates)[indexExcludeVar],"(negative Fisher Information )")

    estimates <- estimates[-indexExcludeVar]

  }

  if (exists ("fisher_info")){  # if dim == 0 no CI is calculated

    standerror <- sqrt (diag (fisher_info))
    numberVarCI <- length (estimates)  # used as an index
    lenghtAlpha <- length (alpha) # used as an index

    # Calculation of CIs  ---------------------------

    # Alpha can be a vector of different significance values
    upper <- matrix (c(rep (NA, numberVarCI * lenghtAlpha)), ncol = lenghtAlpha)
    lower <- matrix (c(rep (NA, numberVarCI * lenghtAlpha)), ncol = lenghtAlpha)

    colnames(upper) <- c(1 - alpha[c(1:lenghtAlpha)] / 2)
    colnames(lower) <- c(alpha[c(1:lenghtAlpha)] / 2)

    # Computation of CIs. The sigmas are chi-squared distributed
    for (i in 1:numberVarCI && dim (fisher_info)[1] > 0){
      if( any (names (estimates[1]) == c("sigma_u", "sigma_v"))){
        upper[i, ] <- df * estimates[i] / qchisq (alpha / 2, df)
        lower[i, ] <- df * estimates[i] / qchisq (1 - alpha / 2, df)
      } else {
        upper[i, ] <- estimates[i] + qnorm (1 - alpha / 2)  * standerror[i]
        lower[i, ] <- estimates[i] + qnorm (alpha / 2) * standerror[i]
      }
    }

    # Create output data.frame for CIs  ---------------------------

    interval <- data.frame (value = estimates, lower = lower, upper = upper, standerror = standerror)
    return (interval)

  } else { # We don`t compute any CIs
    # User is informed that the Hessian Matrix is not valid.
    # It does not stop() terminate other functions applying this function.
    warning ("Hessian Matrix is singular / indefinite. Could not calculate CIs")
  }
}
