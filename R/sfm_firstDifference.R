#' < First-difference transformation & calculation of log.likelihood >
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @param xv is a n*t x k matrix (explantatory variables)
#' @param z is a n*t x r matrix (inefficency determinants)
#' @param y is a n*t x 1 vector (response)
#' @param N is a integer (n - panels)
#' @param Time is a integer (observations per panel)
#' @param group an optional vector specifying the panels to be used in the fitting process.
#' @param mu is a integer (mean of the truncated normal distribution of the inefficency)
#' @param optim is a boolean (set F to obtain a list of model variables.
#'     T to obtain the -sum of log.likelihood)
#' @return If optim = T the log.likelihood is returned of all panels.
#'     If optim = F the model fit is returned including all important model variables.

SFM.firstDiff <- function(par = c(sigma_u, sigma_v, beta = c(), delta = c()),
                          xv, y, z, N = NULL,  Time = NULL, group = NULL, mu=0, optim = F){

  K <- dim (as.matrix (xv))[2]  # K beta variables
  R <- dim (as.matrix (z))[2]  # R delta variables

  if(length (Time) == 1){
    Time <- rep (Time, N)
  }

  cumTime <- c(0, cumsum(Time)) # used for the index of the different panels

  # First Difference of x
  x.diff <- matrix(c(rep(NA, sum(Time - 1) * K)), ncol = K)  # -1 for first diff

  for(u in 1:K){
    difference <- c()
    for(i in 1:N){
      difference <- c(difference, diff(xv[(cumTime[i] + 1) : cumTime[i + 1], u]))
    }
    x.diff[, u] <- difference
  }

  # First Difference of z
  # z.diff <- matrix(c(rep(NA, sum(Time - 1) * R)), ncol = R)  # -1 for first diff
  #
  # for(u in 1:R){
  #   difference <- c()
  #   for(i in 1:N){
  #     difference <- c(difference, diff(z[(cumTime[i] + 1) : cumTime[i + 1], u]))
  #   }
  #   z.diff[, u] <- difference
  # }

  # First DIfference of y
  y.diff <- c()
  for(i in 1:N){
    y.diff <- c(y.diff, diff(y[(cumTime[i]+1):cumTime[i+1], ]))
  }

  epsilon <- y.diff - x.diff %*% par[3:(3+K-1)]  # is a (N * (Time-1)) x 1 vector

  # TODO() i have no idea if this is correct. it should, because the dimensions fit...
  h <- exp (as.matrix (z) %*% par[(4+K-1):(4+K+R-2)])  # R-delta coefficients are used
  # TODO(Clemens): Check if this is correct.

  # splits the vector to N-lists of Time-1 observations
  cumTimeDiff <- c(0, cumsum(Time-1))
  h.diff <- lapply (unname (split (h, findInterval (seq_along (h), cumTime, left.open = TRUE))), diff)
  eps.diff <- unname (split (epsilon, findInterval (seq_along (epsilon), cumTimeDiff, left.open = TRUE)))
  # TODO(Clemens): lapply(..., diff)  not required?

  # if all Time entries are equal we can save time in the next for loop
  # if (diff(range(Time)) < .Machine$double.eps ^ 0.5){
  #   C <- diag(sqrt(par[2]), Time[1])
  #   D <- diff(C)
  #   PI <- D %*% t(D)
  #   invPI <- solve(PI)
  # }

  log.ll      <- c(rep (NA, N))
  mu_2star    <- c(rep (NA, N))
  sigma_2star <- c(rep (NA, N))

  for (i in 1:N){

    # TODO(Clemens) This is inefficient.
    C <- diag(sqrt(par[2]), Time[i])
    D <- diff(C)
    PI <- D %*% t(D)
    invPI <- solve(PI)

    mu_2star[i] <- (mu/par[1] - t(eps.diff[[i]]) %*% invPI %*% h.diff[[i]]) /
      (t(h.diff[[i]]) %*% invPI %*% h.diff[[i]] + 1 / par[1])

    sigma_2star[i] <- 1 / ((t(h.diff[[i]]) %*% invPI %*% h.diff[[i]] + 1 / par[1]))

    log.ll[i] <-  -0.5 * (Time[i] - 1) * log(2 * pi) - 0.5 * log(Time[i]) - 0.5 * (Time[i] - 1) * log(par[2]) -
      0.5 * t(eps.diff[[i]]) %*% invPI %*% eps.diff[[i]] +
      0.5 * ((mu_2star[i]^2) / sigma_2star[i] - (mu^2) / par[1]) +
      log (sqrt (sigma_2star[i]) * pnorm(mu_2star[i] /  sqrt (sigma_2star[i]))) -
      log (sqrt (par[1]) * pnorm (mu / sqrt (par[1])))
  }  # TODO(Clemens) -> check again ll

  # Return values ---------------------------
  if (optim == T){
    # If SFM.within() is called by an optimizer, we need a negative sum of log.ll
    return (sum ((log.ll)*-1))

  } else {
    ret.list         <- list(x.diff = x.diff, y.diff, PI, eps.diff, h.diff, h,
                             mu_2star, sigma_2star, log.ll)
    names (ret.list) <- c("x.wthn","y.wthn", "PI", "eps.wthn", "h", "h.wthn",
                          "mu_2star", "sigma_2star", "log.ll" )
    return (ret.list)
  }

}

