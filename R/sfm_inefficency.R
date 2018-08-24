#' < calculate the inefficiency of each individuum >
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < returns the inefficiency for each individuum >
#' @examples
#' < an example >
#' @export

SFM.inindex <- function(h, sigma2star, mu2star, N, Time, method, cumTime){

  sigmaStar <- sqrt (sigma2star)
  h <- as.matrix (h)

  # calculates the mean inefficency per panel
  inefficencyIndex <- lapply (1:N, function(x) mean (h[(cumTime[x] + 1) : cumTime[x + 1]] *
                                               (mu2star[x] + (dnorm (mu2star[x] / sigmaStar[x]) *
                                               sigmaStar[x]) / pnorm(mu2star[x]/ sigmaStar[x]))))

  inefficencyIndex <- as.matrix (unlist (inefficencyIndex))
  # TODO(): add exp inefficency and maybe plot? for evaluation of quality of results?
  return(inefficencyIndex)
}



