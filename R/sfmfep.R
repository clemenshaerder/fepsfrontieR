#' sfmfep Estimates a Stochastic Frontier Model for Fixed-Effects using
#' within transformation
#' @param formula an object of class "fomrula" in the form of
#' y ~ x1 + ... + x_k + (z1 + ... + z_r). The details of model specification are given under Details
#' @param data an optional data frame, list or environment (or object coercible by
#' as.data.frame to a data frame) containing the variables in the model.
#' @param group an optional vector specifying the panels to be used in the fitting process.
#' @param method a required string specifying the method ("within" or "firstdiff").
#' @param N an optional integer specifying the total amount of panels in the data set.
#' @param Time an optional integer specifying the amount of observations per panel.
#' @param mu is the mean of a truncated normal distribution of the stochastic inefficencys.
#' @param Bootstrap
#' @param sigmaCI is an optional vector specifying the significance values of the confidence intervals
#' for the MLE estimates (based no the Hessian).
#' @param estimate T or F specifies if "myPar" is used as starting point of the estimation,
#' or if "myPar" is used to calculate a given stochastic frontier model.
#' @param myPar is a vecor which has to be entered in the following order:
#' c(sigma_v, sigma_u, beta = c(), delta = c())
#'
#' @return
#' A S3 object is returned including: ....
#' @examples
#' < an example where a data set is created via sfm.generate>
#' @export

sfmfep <- function(formula, data, group = NULL, N = NULL, Time = NULL,
                   mu = 0,  sigmaCI = 0.05, estimate = T, method = "firstdiff", bootstrap = F, B = NULL,
                   myPar = c(sigma_u = NULL, sigma_v = NULL, beta = c(NULL), delta = c(NULL))){


  # Error handling of input data & formula  ---------------------------

  call <- match.call ()

  # Select the data from the "data" input according to applied "formula"
  if ((is.data.frame (data) || is.matrix (data)) == F){
    stop ("Your data should be a data.frame or matrix")
  }

  if (!any (method == c("within", "firstdiff"))){
    stop("A method correct method must be chosen. within or firstdiff")
  }

  # Tests if a correct formula has been plugged in (i.e. y ~ x...)
  try (vars.ex <- all.vars (formula), silent = T)
  if (exists ("vars.ex") == F){  # if the "wiggle ~" is missing, vars.ex will be of lenght 0
    stop ("Formula inadequate.
          Use following notation y ~ x.1 + ... + x.k + (z1 + ... +z.r)")
  }

  # Tests if all variables in the formula match with a column name
  try (sel.data <- data[vars.ex], silent = T)
  if (exists ("sel.data") == F){
    stop ("Couldnt match input *formula* with colnames of input *data*.
          Use the same colnames as in your input file.")
  }

  # Tests if any of the model variables are not numeric
  if (any (sapply (sel.data, function (sel.data) !is.numeric (sel.data)))){
    stop ("All variables must be numeric.")
  }

  # Tests if any NAs exist in the data sel.data
  dataInvalid <- sapply (data, function (sel.data) sum (is.na (sel.data)) +
                                                   sum (is.infinite (sel.data)) +
                                                   sum (is.null (sel.data)))

  if (sum (dataInvalid) > 0){
    print ("Function Stopped. You have NAs/Inf/NULL in your selected data. Summary:")
    stop (print (dataInvalid))
  }

  # Tests if myPar is correctly entered if it is not NULL
  if (!is.null(myPar)){
    if (length (myPar) != (2 + K + R)){  # 2 sigmas + k betas & r deltas
      stop ("Could not perform estimation.
                  A starting point for the estimation must be provided for every parameter.")
    }
    if (!is.numeric (myPar) || myPar[1] <= 0 || myPar[2] <= 0){
      stop ("Could not perform estimation.
                  The starting points must be numeric &| sigmas must be <= 0")
    }
  }

  # Tests if B is correctly specified if bootstrapping is performed
  if (bootstrap == T){
    if (B <= 0 || !is.numeric(B) || length(B) > 1){
      stop("Could not perform calculations.
           B must be an integer > 0")
    }
  }

  # Data Assignment: Translate formula, get N & Time and assign to y, x, z  ---------------------------

  # Get amount of z variables (in curly brackets)
  formula <- as.character (formula)
  zForm <- stringr::str_extract (formula, "(?<=\\().*?(?=\\))")[3]  # extracts everything in curly brackets.
  extractZ <- as.formula (paste (formula[2], formula[1], zForm))

  R <- length (all.vars (extractZ))-1  # -1, as extractZ is of form y ~ z1 + ... + zr
  if (R <= 0){
    stop ("No required z variable were defined in the formula.
          Please add at least one z variable in your formula in brackets (z_i)")
  }
  totCountVar <- dim (sel.data)[2]  # total amount of variables
  K <- totCountVar - R - 1 # all variables - r Z-variables - (1) y-variable = K x variables

  # N & T, or group must be assign, else we can not compute N & T
  # First, we check if no option is defined

  if ((is.null (N) | is.null (Time)) & is.null (group) ){
    stop ("You have to either specify N = panels & Time = obs. per panel
          or provide a group column")
  # Second, we check if N & T are both well definded and if group is not
  } else if ((is.double (N) && is.double (Time)) && !is.character(group)){

    N.input <- N; Time.input <- Time

    if ( sum (N.input * Time.input) != dim (sel.data)[1]) {
      stop ("Your input (N & T) doesn`t match the data dimensions.
            Calculations could not be performed.")
    }

  } else {  # Third if group is specified, we check if this group exists as column name
    if (is.null (group) || try (exists (group, data) == F, silent = T)){
      stop ("Couldnt match input *group* with colnames")
    } else {  # Finally, if either only group is chosen or all three we use group to get N & T
      sel.data <- cbind (sel.data, data[group])  # adds group to the selected data
      sel.data <- sel.data %>% arrange_(.dots = group)  # sorts by group
      N.input <- dim (table (sel.data[group]))[1]
      Time.input <- table (sel.data[group])
    }
  }

  y.dat <- as.matrix(sel.data[, 1])  # y is always the first value in the formula
  x.dat <- as.matrix(sel.data[, 2:(1+K)])
  z.dat <- as.matrix(sel.data[, (1+K+1):(1+K+R)])


  # Estimation  ---------------------------

  #  lower boundary for optimization
  l.int <- c(0.0001, 0.0001, rep(-Inf, K), rep(-Inf, R))  # Variation can not be negative

  # TODO() - hier ist closure fehler
  if (estimate == T){
    if (is.null(myPar) == T){  # we generate appropriate starting values

      beta.start  <- solve (t(x.dat) %*% x.dat) %*% t(x.dat) %*% y.dat  # OLS - within model for beta
      delta.start <- solve (t(z.dat) %*% z.dat) %*% t(z.dat) %*% y.dat

      e <- y.dat - x.dat %*% beta.start
      startSigma <- (t(e) %*% e) / (N.input * min (Time.input) - (K+R))  # OLS for both sigmas

      myPar     <- c(sigma_u = startSigma,
                     sigma_v = startSigma,
                     beta = beta.start,
                     delta = delta.start)
      if (bootstrap == F) {
        if (method == "within"){
          optim.SFM <- nlminb (objective = SFM.within,
                               lower = l.int,
                               hessian = T,
                               start = myPar,
                               Time = Time.input,
                               N = N.input,
                               xv = x.dat, y = y.dat, z = z.dat,
                               mu = mu,
                               optim = T)
        } else {
          optim.SFM <- nlminb (objective = SFM.firstDiff,
                               lower = l.int,
                               hessian = T,
                               start = myPar,
                               Time = Time.input,
                               N = N.input,
                               xv = x.dat, y = y.dat, z = z.dat,
                               mu = mu,
                               optim = T)
        }
      } else {
        optim.SFM <- SFM.bootstrap(y = y.dat, xv = x.dat, z = z.dat,
                                   mu = mu,
                                   N = N.input,
                                   Time = Time.input,
                                   myPar = myPar,
                                   R = R, K = K,
                                   B = B,
                                   method = method,
                                   lowerInt = l.int,
                                   sigmaCI = sigmaCI)
      }
    } else {
      if (bootstrap == F){
        if (method == "within"){
          optim.SFM <- nlminb(objective = SFM.within,
                              start = c(myPar),
                              lower = l.int,
                              mu = mu,
                              Time = Time.input,
                              N = N.input,
                              xv = x.dat, y = y.dat, z = z.dat,
                              optim = T)
        } else {  # use first-difference method
          optim.SFM <- nlminb(objective = SFM.firstDiff,
                              start = c(myPar),
                              lower = l.int,
                              mu = mu,
                              Time = Time.input,
                              N = N.input,
                              xv = x.dat, y = y.dat, z = z.dat,
                              optim = T)
        }
      } else {
        optim.SFM <- SFM.bootstrap(y = y.dat, xv = x.dat, z = z.dat,
                                   mu = mu,
                                   N = N.input,
                                   Time = Time.input,
                                   myPar = myPar,
                                   method = method,
                                   R = R, K = K,
                                   B = B,
                                   lowerInt = l.int,
                                   sigmaCI = sigmaCI)
      }
    }
  } else { # we dont estimate (estimate = F)
      # TODO(Clemens): Maybe extend this for Bootstrapping too.
      optim.SFM <- NULL
      optim.SFM$par <- myPar  # provided parameters are used as parameters
  }

  # derive the Hessian Matrix based on estimates from the optimization
  if (bootstrap == F){
    hes <- numDeriv::hessian(SFM.within,
                           method = "Richardson",
                           x = optim.SFM$par,
                           xv = x.dat, y = y.dat, z = z.dat,
                           N = N.input,
                           Time = Time.input,
                           mu = mu,
                           optim = T)  # Note optim = T

    #standard errors for estimates
    standerror <- sqrt (diag (solve (hes))) # TODO() Check confidence Interval -> we do that twice (and only 1 time with error handler)

  }



  # Fit the model based on the estimation  ---------------------------

  # TODO() Do we actually need that???
  estimate.sigma_u  <- optim.SFM$par[1]
  estimate.sigma_v  <- optim.SFM$par[2]
  estimate.beta     <- optim.SFM$par[3:(3+K-1)]
  estimate.delta    <- optim.SFM$par[(3+K):(3+R-1)]

  # Fit the specified model with the estimates
  if (method == "within"){
    ret.list <- SFM.within(optim = F,  # Note optim = F
                           N = N.input,
                           Time = Time.input,
                           xv = x.dat , y = y.dat ,z = z.dat,
                           par = optim.SFM$par)
  } else {
    ret.list <- SFM.firstDiff(optim = F,  # Note optim = F
                              N = N.input,
                              Time = Time.input,
                              xv = x.dat , y = y.dat ,z = z.dat,
                              par = optim.SFM$par)
  }

  if (estimate == F || bootstrap == T){  # log.ll is obtained from the model as nlminb is not primarilyused.
    optim.SFM$objective <- sum(ret.list$log.ll*-1)
  }


  # Calculated Confidence Intervals  ---------------------------

  # sigmaCI can be a vector of significance levels
  # TODO(Clemens): only calculated if bootstrapping is not performed
  if (bootstrap == F){
    if ( (!any (sigmaCI <= 0 | sigmaCI > 1)) && !is.null (sigmaCI) && !is.nan (sigmaCI) ){
      # A data frame is returned
      conf.Interval <- SFM.CI(estimates = optim.SFM$par, hessianMatrix = hes, alpha = sigmaCI)
    } else if (is.null (sigmaCI)){

    } else if (any (sigmaCI <= 0 | sigmaCI > 1) || is.nan (sigmaCI)){
      cat ("Could not compute Confidence Intervals due to invalid input (sigmaCI must be between [0, 1]")
    }
  }

  # Inefficency Index for each Panel  ---------------------------

  # TODO(Oli): extend inefficency to Time as an vector (check sfm_within for that)
  # inefficency <- SFM.inindex(h = ret.list$h,  # Note h is not within transformed
  #                            sigma2star = ret.list$sigma_2star,
  #                            mu2star = ret.list$mu_2star,
  #                            N = N.input,
  #                            Time = Time.input)

  # Recover Fixed Effects (alpha) for each Panel  ---------------------------

  # TODO(Oli): extend alpha to Time as an vector (check sfm_within for that)
  # alpha <- SFM.alpha(beta = estimate.beta,
  #                    mu = mu,
  #                    sigma_u = estimate.sigma_u,
  #                    sigma_v = estimate.sigma_v,
  #                    h = ret.list$h,
  #                    x = x.dat,
  #                    y = y.dat,
  #                    epsilon = ret.list$eps.wthn,
  #                    N = N.input,
  #                    Time = Time.input)

  # Model Selection Criterion  ---------------------------

  # Adjusted formula with "* -1" as the used objective is the "-"log.ll (needed for optimization)
  AIC <- -2 * -1 * optim.SFM$objective + 2 * length (optim.SFM$par)
  BIC <- -2 * -1 * optim.SFM$objective + length (optim.SFM$par) * dim (y.dat)[1]


  # Output  ---------------------------

  # TODO(Oli): add option if some values are not created like conf.Interval
  # -> we have currently one fail in unit tester due   to that
  res <- list (call = call, par = myPar,
               estimates = optim.SFM$par , AIC = AIC, BIC = BIC, estimate = estimate,
              # ci = conf.Interval
              # ret.list = ret.list,
              # alpha = alpha
              # inefficency = inefficency
              # hessian = hes
              # standerror = standerror
              contrasts = c(attr (myPar, "names")),
              likeihood= optim.SFM$objective)

  class (res) <- c(res$class, "sfmfep")
  res
}
