#' @title Within
#' @description Performs a within transformation to a Stochastic Frontier Model.
#' By within-transformation, the sample mean of each panel is subtracted
#' from every observation in the panel. The transformation thus removes
#' the time-invariant individual effect from the model
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
#' @return If optim = T the log.likelihood is returned of all panels.
#'     If optim = F the model fit is returned including all important model variables.
#' @export


SFM.within <- function(par = c(sigma_u, sigma_v, beta = c(), delta = c()),
                       cumTime, xv, y, z, N = NULL,  Time = NULL, mu=0,
                       optim = F, K = NULL, R = NULL, seqN = 1:N){

  # Within Transformations ---------------------------

  # compute interval for within transformations of x, epsilon & h
  splitInterval <- findInterval (seq_along (y[, ]), cumTime, left.open = TRUE)


  # Within transformation of X
  # A vector of length Ti (obseravations per panel) is created and binded
  # to the existing vector until the vector is of length N*(Ti).
  # We do this for K parameters.
  x.wthn <- lapply(1:K, function(k)
    by (xv[, k], INDICES = splitInterval, FUN = function(x) x - mean(x)
    ))
  x.wthn <- matrix (unname (unlist (x.wthn)), ncol = K)

  # Within transformation of Y
  # Same transformation procedure as for x
  repYMeans <- c()
  repYMeans <- lapply (seqN, function(x)
    repYMeans <-
      c(repYMeans, rep (mean (y[(cumTime[x] + 1):cumTime[x + 1],]), Time[x])))
  repYMeans <- matrix (unlist (repYMeans))

  y.wthn <- y - repYMeans

  # Computes the residuals of the centered response & explanatory variables
  # based on beta-estimates
  epsilon <- as.matrix (y.wthn - x.wthn %*% par[3:(3+K-1)])

  # Within transformation of epsilon is not required (already mean 0).
  # Note, that epsilon is a vector that needs to be split to lists
  # for the likelihood computation
  eps.wthn <- unname (split (epsilon, splitInterval))

  # Within Transformation of h
  # An exponential function is applied on the z inefficency determinants
  h      <- exp (as.matrix (z) %*% par[(4+K-1):(4+K+R-2)])  # R-delta coefficients are used
  h.wthn <- lapply (unname (split (h, splitInterval)),
                     function(x) x - rep (mean (x), length (x)))


  # Log-Likelihood computation for each panel ---------------------------

  # PI-Matrix is the Variance-Covariance Matrix of v
  # It is a Time_i x Time_i symmetric matrix.
  #Thus we can just compute it once for max Time
  mTime <- max (Time)

  PI <-
    par[2] * (diag (mTime) - 1 / mTime * (matrix (c(rep (1, mTime * mTime
    )), ncol = mTime)))

  try(gPI <- ginv(PI), silent = T)
  if (!exists("gPI")) {
    stop (
      "Could not calculate log.likelihood.
      SVD of the g-Inverse of PI failed.
      PI is sigma_v * M (M is the TxT orthogonal projection matrix).
      Try instead *method = firstdiff*."
    )
  }

  # PI matrix is adjusted for each panel to Ti x Ti
  itPI <- lapply(Time, function(x) gPI[(1:x), (1:x)])

  # calculates "mu two star" needed for the log-likelihood for each panel
  mu_2star <- lapply (seqN, function(x)
    (mu / par[1] - t(eps.wthn[[x]]) %*% itPI[[x]] %*% h.wthn[[x]]) /
      (t(h.wthn[[x]]) %*% itPI[[x]] %*% h.wthn[[x]] + 1 / par[1]))
  mu_2star <- matrix (unlist (mu_2star))

  # calculates "sigma two star" needed for the log-likelihood for each panel
  sigma_2star <- lapply (seqN, function(x)
      1 / ((t(h.wthn[[x]]) %*% itPI[[x]] %*% h.wthn[[x]] + 1 / par[1])))
  sigma_2star <- matrix (unlist (sigma_2star))

  # log-likelihood for N panels
  log.ll <- lapply (seqN, function(x)
    - 0.5 * (Time[x] - 1) * log(2 * pi) -
      0.5 * (Time[x] - 1) * log(par[2]) -
      0.5 * t(eps.wthn[[x]]) %*% itPI[[x]] %*% eps.wthn[[x]] +
      0.5 * ((mu_2star[x] ^ 2) / sigma_2star[x] - (mu ^ 2) / par[1]) +
             log (sqrt (sigma_2star[x]) * pnorm(mu_2star[x] /  sqrt (sigma_2star[x]))) -
      log (sqrt (par[1]) * pnorm (mu / sqrt (par[1]))))
  log.ll <- matrix (unlist (log.ll))


  # Return values ---------------------------
  # "else" for cases where optimization is not desired and another optimization
  # procedure may be used

  if (optim == T) {
    # If SFM.within() is called by an optimizer, we need a negative sum of log.ll
    return (sum ((log.ll) * -1))
  } else {
    ret.list <-
      list (
        x.trans = x.wthn,
        y.trans = y.wthn,
        PI,
        eps.trans = eps.wthn,
        h.trans = h.wthn,
        h = h,
        mu_2star = mu_2star,
        sigma_2star = sigma_2star,
        log.ll = log.ll
      )
    names (ret.list) <-
      c(
        "x.trans",
        "y.trans",
        "PI",
        "eps.trans",
        "h.trans",
        "h",
        "mu_2star",
        "sigma_2star",
        "log.ll"
        )
    return (ret.list)
  }
}
