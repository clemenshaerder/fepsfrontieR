#' < calculate the inefficiency of each individuum >
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @return < returns the inefficiency for each individuum >
#' @examples
#' < an example >
#' @export

SFM.inindex.unbalanced <- function(h, sigma2star, mu2star, N, Time, method){

  # TODO(Oli) add panelname to each ineff.index
  if(length (Time) == 1){
    Time <- rep (Time, N)
  }

  sigmaStar <- sqrt (sigma2star)
  cumTime <- c(0, cumsum (Time))
  h <- as.matrix(h)
  inIndex_i <- c()
  inefficencyIndex <- c()

  for (i in 1:N){
    inIndex_i <- c(inIndex_i, h[(cumTime[i] + 1) : cumTime[i + 1]] *
                  (mu2star[i] + (dnorm(mu2star[i] / sigmaStar[i]) *
                                  sigmaStar[i]) / pnorm(mu2star[i]/ sigmaStar[i])))
  }

  for(j in 1:N){
    inefficencyIndex <- c(inefficencyIndex, mean( inIndex_i[(cumTime[j] + 1) : cumTime[j + 1]]))
  }

  inefficencyIndex <- as.matrix(inefficencyIndex)

  return(inefficencyIndex)
}



