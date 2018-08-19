#' Performs a within transformation to a Stochastic Frontier Model
#'
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @param xv is a n*t x k matrix (explantatory variables)
#' @param z is a n*t x r matrix (inefficency determinants)
#' @param y is a n*t x 1 vector (response)
#' @param N is a integer (n - panels)
#' @param Time is a integer (observations per panel)
#' @param mu is a number (mean of the truncated normal distribution of the inefficency)
#' @param optim is a boolean (set F to obtain a list of model variables.
#'     T to obtain the -sum of log.likelihood)
#' @param group an optional vector specifying the panels to be used in the fitting process.
#' @return If optim = T the log.likelihood is returned of all panels.
#'     If optim = F the model fit is returned including all important model variables.
#' @export


SFM.within <- function(par = c(sigma_u, sigma_v, beta = c(), delta = c()),
                       xv, y, z, N = NULL,  Time = NULL, group = NULL, mu=0, optim = F){

  K <- dim (as.matrix (xv))[2]  # K beta variables
  R <- dim (as.matrix (z))[2]  # R delta variables

  # Within Transformations ---------------------------

  # In case of a balanced panel each N has Time observations
  if(length (Time) == 1){
    Time <- rep (Time, N)
  }

  cumTime <- c(0, cumsum (Time)) # used for the index of the variables

  # Within transformation of X
  x.wthn <- matrix(c(rep (NA, sum (Time) * K)), ncol = K)
  for(u in 1:K){  # for k explenatory variables
    repMeans <- c()
    for(i in 1:N){  # do it for each panel
      repMeans <- c(repMeans, rep(mean (xv[(cumTime[i] + 1) : cumTime[i + 1],u]), Time[i]))
    }
    x.wthn[, u] <- (xv[, u] - repMeans)
  }
  x.wthn <- as.matrix (x.wthn, ncol = K)

  # Within transformation of Y
  repYMeans <- c()
  for(i in 1:N){
    repYMeans <- c(repYMeans, rep(mean (y[(cumTime[i]+1):cumTime[i+1], ] ), Time[i]))
  }
  y.wthn <- as.matrix(y - repYMeans)


  # Computes the residuals of the centered response & explanatory variables based on beta-estimates
  epsilon <- y.wthn - x.wthn %*% par[3:(3+K-1)]
  splitInterval <- findInterval (seq_along (epsilon), cumTime, left.open = TRUE)

  # Within transformation of epsilon is not required (already mean 0).
  # Note, that epsilon is a vector that needs to be split to lists for the likelihood computation

  eps.wthn <- unname (split (epsilon, splitInterval))

  # Within Transformation of h
  # An exponential function is applied on the z inefficency determinants
  h      <- exp (as.matrix (z) %*% par[(4+K-1):(4+K+R-2)])  # R-delta coefficients are used
  h.wthn <- lapply (unname (split (h, splitInterval)),
                     function(x) x - rep (mean (x), length (x)))

  # Log-Likelihood computation for each panel ---------------------------

  # PI-Matrix is the Variance-Covariance Matrix of v
  # It is a Time_i x Time_i symmetric matrix. Thus we can just compute it once for max Time
  mTime <- max(Time)

  PI <- par[2] * (diag (mTime) - 1/mTime * (matrix (c(rep (1, mTime * mTime)), ncol = mTime)))
  try(gPI <- MASS::ginv(PI), silent=T)
  if (!exists("gPI")){
    stop ("Could not calculate log.likelihood.
          SVD of the g-Inverse of PI failed.
          PI is sigma_v * M (M is the TxT orthogonal projection matrix).")
  }

  log.ll      <- c(rep (NA, N))
  mu_2star    <- c(rep (NA, N))
  sigma_2star <- c(rep (NA, N))

  # the likelihood for each panel is computed applying the list entries of e.wthn & h.wthn.
  for (i in 1:N){

    # for each panel we adapt the dimensions of the gPI which was calculated based on the max(Time)
    itPI <- gPI[(1:Time[i]), (1:Time[i])]

    mu_2star[i] <- (mu/par[1] - t(eps.wthn[[i]]) %*% itPI %*% h.wthn[[i]]) /
                   (t(h.wthn[[i]]) %*% itPI %*% h.wthn[[i]] + 1 / par[1])

    sigma_2star[i] <- 1 / ((t(h.wthn[[i]]) %*% itPI %*% h.wthn[[i]] + 1 / par[1]))

    log.ll[i] <-  -0.5 * (Time[i] - 1) * log(2 * pi) - 0.5 * (Time[i] - 1) * log(par[2]) -
                  0.5 * t(eps.wthn[[i]]) %*% itPI %*% eps.wthn[[i]] +
                  0.5 * ((mu_2star[i]^2) / sigma_2star[i] - (mu^2) / par[1]) +
                  log (sqrt (sigma_2star[i]) * pnorm(mu_2star[i] /  sqrt (sigma_2star[i]))) -
                  log (sqrt (par[1]) * pnorm (mu / sqrt (par[1])))
  }

  # Return values ---------------------------
  if (optim == T){
    # If SFM.within() is called by an optimizer, we need a negative sum of log.ll
    return (sum ((log.ll)*-1))

  } else {
    ret.list         <- list(x.wthn = x.wthn, y.wthn, PI, eps.wthn, h.wthn, h,
                             mu_2star, sigma_2star, log.ll)
    names (ret.list) <- c("x.wthn","y.wthn", "PI", "eps.wthn", "h", "h.wthn",
                          "mu_2star", "sigma_2star", "log.ll" )
    return (ret.list)
  }

}
