#' @title First Differnces Transformation of a Stochastic Frontier Model
#' @description Performs a first-difference transformation to a Stochastic Frontier Model
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @param cumTime ia a vector of the cumulated times of the Time vector.
#'     It serves as an index for computation.
#' @param xv is a n*t x k matrix (explantatory variables)
#' @param y is a n*t x 1 vector (response)
#' @param z is a n*t x r matrix (inefficency determinants)
#' @param N is a integer (n - panels)
#' @param Time is a integer (observations per panel)
#' @param mu is a integer (mean of the truncated normal distribution of the inefficency)
#' @param optim is a boolean (set F to obtain a list of model variables.
#'     T to obtain the -sum of log.likelihood)
#' @param K is an integer (# of xv variables)
#' @param R is an integer (# of z variables)
#' @param seqN is a sequence from 1 to N
#' @return If optim = T the log.likelihood is returned of all panels.
#'     If optim = F the model fit is returned including all important model variables.

SFM.firstDiff <- function(par = c(sigma_u, sigma_v, beta = c(), delta = c()),
                          cumTime, xv, y, z, N = NULL,  Time = NULL, mu=0,
                          optim = F, K = NULL, R = NULL, seqN = 1:N){


  # First-difference Transformations ---------------------------

  # y has always the same lenght as h & x, thus we can reduce calculation
  # time by calculating the index sequence just once in the function
  seqDifference <- findInterval ( seq_along (y), cumTime, left.open = TRUE)

  # First-Difference of X
  # For every col of x the first difference of each panel is calculated,
  # resulting in a N * (Ti - 1) x K matrix
  x.diff <- matrix(c(rep (NA, sum (Time - 1) * K)), ncol = K)  # -1 for first diff
  diff <- replicate (K, list (), simplify = F)  # create R bootstrap samples

  diff <- lapply(1:K, function(k) lapply ( unname (split (xv[, k], seqDifference)), diff))
  x.diff <- matrix (unlist (diff), ncol=K)

  # First-Difference of y
  y.diff <- lapply ( unname (split (y, seqDifference)), diff)
  y.diff <- matrix (unlist (y.diff, recursive = T), ncol = 1)

  h <- exp (as.matrix (z) %*% par[(4+K-1):(4+K+R-2)])  # R-delta coefficients are used

  # First difference of h
  h.diff <- lapply (unname (split (h, seqDifference)), diff)

  # Computes the residuals of the centered response & explanatory variables
  # based on beta-estimates. Epsilon is already transformed.
  epsilon <- y.diff - x.diff %*% par[3:(3+K-1)]  # is a (N * (Time-1)) x 1 vector

  cumTimeDiff <- c(0, cumsum(Time-1))
  # create lists to continue computations. epsilon need a different sequence.
  eps.diff <-
    unname (split ( epsilon,
      findInterval (seq_along (epsilon), cumTimeDiff, left.open = TRUE)
    ))

  # notation of the covariance matrix of "v" is labeled here and following as PI
  # Wang & Ho uses the notation sigma
  # Computation of the Inverse of the PI matrix for each panel -----------------

  # if all Time entries are equal, we can calculate one invPI to
  # significantly increase the efficency
  if (diff(range(Time)) < .Machine$double.eps ^ 0.5){
    C     <- diag(sqrt(par[2]), Time[1])
    D     <- diff(C)
    PI    <- D %*% t(D)
    invPI <- solve(PI)
    invPI <- replicate(N, invPI, simplify = F)  # create N PI-inverse as a list.

  } else {  # Time is unequal and thus we have unequal panels with differen PI-inverse
    D     <- lapply(Time, function(x) diff(diag(sqrt(par[2]), x)))
    invPI <- lapply(seqN, function(x) solve(D[[x]] %*% t(D[[x]])))
  }


  # Computation of the log likelihood for each panel ---------------------------

  # calculates "mu two star" needed for the log-likelihood for each panel
  mu_2star <- lapply (seqN, function(x)
    (mu / par[1] - t(eps.diff[[x]]) %*% invPI[[x]] %*% h.diff[[x]]) /
      (t(h.diff[[x]]) %*% invPI[[x]] %*% h.diff[[x]] + 1 / par[1]))
  mu_2star <- matrix (unlist (mu_2star))

  # calculates "sigma two star" needed for the log-likelihood for each panel
  sigma_2star <- lapply(seqN, function(x)
    1 / ((t(h.diff[[x]]) %*% invPI[[x]] %*% h.diff[[x]] + 1 / par[1])))
  sigma_2star <- matrix (unlist (sigma_2star))

  # log-likelihood for N panels
  log.ll <- lapply(seqN, function(x)
    - 0.5 * (Time[x] - 1) * log (2 * pi) -
      0.5 * log (Time[x]) -
      0.5 * (Time[x] - 1) * log (par[2]) -
      0.5 * t(eps.diff[[x]]) %*% invPI[[x]] %*% eps.diff[[x]] +
      0.5 * ((mu_2star[x] ^ 2) / sigma_2star[x] - (mu ^ 2) / par[1]) +
      log (sqrt (sigma_2star[x]) *  pnorm(mu_2star[x] /  sqrt (sigma_2star[x]))) -
      log (sqrt (par[1]) * pnorm (mu / sqrt (par[1]))))

  log.ll <- matrix (unlist (log.ll))


  # Return values ---------------------------
  # "else" for cases where optimization is not desired and another optimization
  # procedure may be used

  if (optim == T) {
    # If SFM.within() is called by an optimizer, we need a negative sum of log.ll
    return (sum ((log.ll) * -1))

  } else {
    ret.list         <- list(
      x.trans = x.diff,
      y.trans = y.diff,
      invPI,
      eps.trans = eps.diff,
      h.trans = h.diff,
      h,
      mu_2star,
      sigma_2star,
      log.ll
    )
    names (ret.list) <-
      c(
        "x.trans",
        "y.trans",
        "invPI",
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

