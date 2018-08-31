#' @title Confidence Intervals
#' @description  Computes Confidence Intervals for the regression coefficients based on the Hessian Matrix
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < Describe what is returned when applying this functinon >
#' @examples
#' < an example >

SFM.CI <- function(estimates, hessianMatrix, alpha, N, Time, df){

  try (fisher_info <- solve (hessianMatrix), silent = T)

  # Exclusion of improper Fisher entries  ---------------------------

  if (exists ("fisher_info") && dim (fisher_info)[1] > 0){

    indexIncludeVar <- which (diag (fisher_info) >= 0)

    if (any (diag (fisher_info) < 0) == T) {
      # If the Hessian Matrix is indefinite, we can not calculate Standard Errors.
      # Occures if eigenvalues of the Hessian are != 0 (estimates are saddle points)
      indexExcludeVar <- which (diag (fisher_info) < 0)
      cat (
        "Could not compute Standard Errors & Confidence Interval for:",
        names(estimates)[indexExcludeVar],
        "(negative Fisher Information )"
      )
    }

    standerror <- rep(NA, length (estimates))

    # Calculation of CIs  ---------------------------

    # Calculates standarderror only for valid optimas
    standerror[indexIncludeVar] <-
      sqrt (diag (fisher_info[indexIncludeVar,
                              indexIncludeVar]))

    # generate matrices for lower & upper bounds.
    upper <-
      matrix (c(rep (NA, length (estimates) * length (alpha))),
              ncol = length (alpha))

    lower <-
      matrix (c(rep (NA, length (estimates) * length (alpha))),
              ncol = length (alpha))

    # colname is the probability bound
    colnames (upper) <- c(1 - alpha[c(1:length (alpha))] / 2)
    colnames (lower) <- c(alpha[c(1:length (alpha))] / 2)

    # Sigmas are asymptotically chi-sq distributed, thus calculated separately
    # "<= 2 " as the first two estimates are always the sigmas
    IndexSigmaCI <- indexIncludeVar[which(indexIncludeVar <= 2)]

    upper[IndexSigmaCI,] <-
      df * estimates[IndexSigmaCI] / qchisq (alpha / 2, df)

    lower[IndexSigmaCI,] <-
      df * estimates[IndexSigmaCI] / qchisq (1 - alpha / 2, df)

    # Calculation of the asymptotic normally distributed betas & deltas
    indexCoeff          <- indexIncludeVar[which(indexIncludeVar > 2)]

    upper[indexCoeff, ] <- estimates[indexCoeff] +
                           qt (1 - alpha / 2, df) * standerror[indexCoeff]

    lower[indexCoeff, ] <- estimates[indexCoeff] +
                           qt (alpha / 2, df) * standerror[indexCoeff]


    # Create output data.frame for CIs  ---------------------------

    interval <- data.frame (
      value = estimates,
      lower = lower,
      upper = upper,
      standerror = standerror
    )
    return (interval)

  } else {
    # else exist(fisher_info) = F -> We don`t compute any CIs
    # User is informed that the Hessian Matrix is not valid.
    # It does not stop() terminate other functions applying this function.
    warning ("Hessian Matrix is singular / indefinite. Could not calculate CIs")
  }
}
