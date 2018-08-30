compare <- function(N, Time, beta, delta,
                    sigma_u, sigma_v, mu, Iteration){

  Time    <- rep(10, N)
  cumTime <- c(0, cumsum (Time))
  K       <- length(beta)
  R       <- length(delta)
  par     <- c(sigma_u,sigma_v,beta,delta)
  lower <- c(0.05,0.05, rep(-Inf,K), rep(-Inf,R))
  upper <-  c(1000,1000,1000,1000,1000,1000)

  estimateW <- matrix(c(rep(NA, (2+K+R)*Iteration)),
                      nrow= 2 + K + R, ncol=Iteration)
  estimateF <- matrix(c(rep(NA, (2+K+R)*Iteration)),
                      nrow= 2 + K + R, ncol=Iteration)

  for(i in 1:Iteration){
    data <- SFM.generate(N = N, Time = min(Time), beta = beta,
                         delta = delta, mu = mu, sigma_u = sigma_u,
                         sigma_v = sigma_v)
    x    <- as.matrix(data[, 1:length(beta)])
    y    <- as.matrix(data[, (length(beta)+1)])
    z    <- as.matrix(data[, (length(beta)+2):(length(beta)+2+length(delta)-1)])

    optim.SFMwithin <- nlminb (
      objective = SFM.within,
      start = c(sigma_u, sigma_v, beta, delta),
      lower = c(0.001, 0.001, rep(-Inf, R + K)),
      mu = mu,
      Time = Time,
      N = N,
      R = length(delta),
      K = length(beta),
      cumTime = cumTime,
      xv = x,
      y = y,
      z = z,
      optim = T
    )
    estimateW[,i] <- optim.SFMwithin$par

    optim.SFMfirstDiff <- nlminb (
      objective = SFM.firstDiff,
      start = c(sigma_u, sigma_v, beta, delta),
      lower = c(0.001, 0.001, rep(-Inf, R + K)),
      mu = mu,
      Time = Time,
      N = N,
      cumTime = cumTime,
      K = K,
      R = R,
      xv = x,
      y = y,
      z = z,
      optim = T
    )
    estimateF[,i] <- optim.SFMfirstDiff$par
  }

  return(cbind(apply(estimateW, 1, sum), apply(estimateF, 1, sum)) / Iteration)
}
# compare this to the estimtes of Wang and HO
compare(N = 100, Time = 10, beta = c(0.5), delta=c(0.5),
        sigma_u = 0.2, sigma_v = 0.1, mu = 0, Iteration = 10)
