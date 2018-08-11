#' < First-difference transformation & calculation of log.likelihood >
#' @param par is a vector of regression coefficients & variance parameters.
#'     1st parameter: sigma_u, 2nd parameter: sigma_v, followed by K beta & R delta coefficients
#' @param xv is a n*t x k matrix (explantatory variables)
#' @param z is a n*t x r matrix (inefficency determinants)
#' @param y is a n*t x 1 vector (response)
#' @param N is a number (n - panels)
#' @param T is a number (observations per panel)
#' @param mu is a number (mean of the truncated normal distribution of the inefficency)
#' @param optim is a boolean (set F to obtain a list of model variables.
#'     T to obtain the -sum of log.likelihood)
#' @return < Describe what is returned when applying this functino >
#' @examples
#' < create an example (look at lm() or something like that) >

SFM.firstDiff <- function(par = c(sigma_u, sigma_v, beta = c(), delta = c()),
                       xv, y, z, N = NULL,  Time = NULL, group = NULL, mu=0, optim = F){

  K <- dim (as.matrix (xv))[2]  # K beta variables
  R <- dim (as.matrix (z))[2]  # R delta variables

  if(length (Time) == 1){
    Time <- rep (Time, N)
  }

  cumTime <- c(0, cumsum(Time)) # used for the index of the variables

  # First Difference of x
  x.diff <- matrix(c(rep(NA, sum(Time - 1) * N)), ncol = K)  # -1 for first diff

  for(u in 1:K){
    difference <- c()
    for(i in 1:N){
      difference <- c(difference, diff(xv[(cumTime[i] + 1) : cumTime[i + 1], u]))
    }
    x.diff[, u] <- difference
  }

  z.diff <- matrix(c(rep(NA, sum(Time - 1) * N)), ncol = R)  # -1 for first diff

  for(u in 1:R){
    difference <- c()
    for(i in 1:N){
      difference <- c(difference, diff(z[(cumTime[i] + 1) : cumTime[i + 1], u]))
    }
    z.diff[, u] <- difference
  }


  y.diff <- c()
  for(i in 1:N){
    y.diff <- c(y.diff, diff(y[(cumTime[i]+1):cumTime[i+1], ]))
  }

  epsilon <- y.diff - x.diff %*% par[3:(3+K-1)]  # is a (N * (Time-1)) x 1 vector

  # TODO() i have no idea if this is correct. it should, because the dimensions fit.. still
  h      <- exp (as.matrix (z.diff) %*% par[(4+K-1):(4+K+R-2)])  # R-delta coefficients are used
  # TODO(Clemens): Check if this is correct.
  cumTimeH <- c(0, cumsum(Time))
  # h.diff <- lapply(unname (split (h, findInterval (seq_along (h), cumTimeH))), diff)
  h.diff <- unname (split (h, findInterval (seq_along (h), cumTimeH)))
  cumTimeE <- c(0, cumsum(Time) - N + 1)
  # TODO(Clemens): lapply(..., diff)  not required?
  eps.diff <- unname (split (epsilon, findInterval (seq_along (epsilon), cumTime)))


  mTime <- max(Time)

  C <- diag(sqrt(par[2]), mTime)
  D <- diff(C)
  PI <- D %*% t(D)

  log.ll      <- c(rep (NA, N))
  mu_2star    <- c(rep (NA, N))
  sigma_2star <- c(rep (NA, N))

  for (i in 1:N){
    itPI <- PI[(1:(Time[i]-1)), (1:(Time[i]-1))]  # for each panel we adapt the dimensions of the gPI
    mu_2star[i] <- (mu/par[1] - t(eps.diff[[i]]) %*% itPI %*% h.diff[[i]]) /
                   (t(h.diff[[i]]) %*% itPI %*% h.diff[[i]] + 1 / par[1])

    sigma_2star[i] <- 1 / ((t(h.diff[[i]]) %*% itPI %*% h.diff[[i]] + 1 / par[1]))

    log.ll[i] <-  -0.5 * (Time[i] - 1) * log(2 * pi) - 0.5 * log(Time[i]) - 0.5 * (Time[i] - 1) * log(par[2]) -
      0.5 * t(eps.diff[[i]]) %*% itPI %*% eps.diff[[i]] +
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

