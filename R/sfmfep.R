#' sfmfep Estimates a Stochastic Frontier Model for Fixed-Effects using
#' within transformation
#' @param formula an object of class "fomrula" in the form of
#' y ~ x1 + ... + x_k + (z1 + ... + z_r). The details of model specification are given under Details
#' @param data an optional data frame, list or environment (or object coercible by
#' as.data.frame to a data frame) containing the variables in the model.
#' @param group an optional vector specifying the panels to be used in the fitting process.
#' @param N an optional integer specifying the total amount of panels in the data set.
#' @param Time an optional integer specifying the amount of observations per panel.
#' @param mu is the mean of a truncated normal distribution of the stochastic inefficencys.
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
                   mu = 0,  sigmaCI = 0.05, estimate = T, method = "within",
                   myPar = c(sigma_u = NULL, sigma_v = NULL, beta = c(NULL), delta = c(NULL))){


  # Error handling of input data & formula  ---------------------------

  call <- match.call()

  # Select the data from the "data" input according to applied "formula"
  if ((is.data.frame(data) || is.matrix(data)) == F){
    stop ("Your data should be a data.frame or matrix")
  }

  # Tests if a correct formula has been plugged in (i.e. y ~ x...)
  try (vars.ex <- all.vars(formula), silent = T)
  if (exists("vars.ex") == F){  # if the "wiggle ~" is missing, vars.ex will be of lenght 0
    stop ("Formula inadequate.
          Use following notation y ~ x.1 + ... + x.k + (z1 + ... +z.r)")
  }

  # Tests if all variables in the formula match with a column name
  try (sel.data <- data[vars.ex], silent = T)
  if (exists("sel.data") == F){
    stop ("Couldnt match input *formula* with colnames of input *data*.
          Use the same colnames as in your input file.")
  }

  # Tests if any of the model variables are not numeric
  if (any(sapply(sel.data, function(sel.data) !is.numeric(sel.data)))){
    stop ("All variables must be numeric.")
  }

  # Tests if any NAs exist in the data sel.data
  dataInvalid <- sapply(data, function(sel.data) sum(is.na(sel.data)) +
                                                 sum(is.infinite(sel.data)) +
                                                 sum(is.null(sel.data)))

  if (sum(dataInvalid) > 0){
    print ("Function Stopped. You have NAs/Inf/NULL in your selected data. Summary:")
    stop (print(dataInvalid))
  }

  # Translate formula, get N & Time and assign to variables  ---------------------------

  # Get amount of z variables (in curly brackets)
  formula <- as.character(formula)
  zForm <- stringr::str_extract(formula, "(?<=\\().*?(?=\\))")[3]  # extracts everything in curly brackets.
  extractZ <- as.formula(paste(formula[2], formula[1], zForm))

  R <- length(all.vars(extractZ))-1  # -1, as extractZ is of form y ~ z1 + ... + zr
  if (R <= 0){
    stop("No required z variable were defined in the formula.
         Please add at least one z variable in your formula in brackets (z_i)")
  }
  totCountVar <- dim(sel.data)[2]  # total amount of variables
  K <- totCountVar - R - 1 # all variables - r Z-variables - (1) y-variable = K x variables

  # N & T, or group must be assign, else we can not compute N & T
  # First, we check if no option is defined

  if ((is.null(N) | is.null(Time)) & is.null(group) ){
    stop ("You have to either specify N = panels & Time = obs. per panel
          or provide a group column")
  # Second, we check if N & T are both well definded and if group is not
  } else if ((is.double (N) && is.double (Time)) && !is.character(group)){

    N.input <- N; Time.input <- Time

    if ( sum(N.input * Time.input) != dim(sel.data)[1]) {
      stop("Your input (N & T) doesn`t match the data dimensions.
           Could not perform calculations.")
    }

  } else {  # Third if group is specified, we check if this group exists as column name
    if (is.null(group) || try(exists(group, data) == F, silent = T)){
      stop ("Couldnt match input *group* with colnames")
    } else {  # Finally, if either only group is chosen or all three we use group to get N & T
      sel.data <- cbind(sel.data, data[group])  # adds group to the selected data
      sel.data <- sel.data %>% arrange_(.dots = group)  # sorts by group
      N.input <- dim(table(sel.data[group]))[1]
      Time.input <- table(sel.data[group])
    }
  }

  y.dat <- as.matrix(sel.data[, 1])  # y is always the first value in the formula
  x.dat <- as.matrix(sel.data[, 2:(1+K)])
  z.dat <- as.matrix(sel.data[, (1+K+1):(1+K+R)])

  # Optimization  ---------------------------
  #  lower boundary for optimization
  l.int <- c(0.0001, 0.0001, rep(-Inf, K), rep(-Inf, R))  # Variation can not be negative

  if (estimate == T){
    if (is.null(myPar) == T){  # we generate appropriate starting values

      beta.start  <- solve (t(x.dat) %*% x.dat) %*% t(x.dat) %*% y.dat  # OLS - within model for beta
      delta.start <- solve (t(z.dat) %*% z.dat) %*% t(z.dat) %*% y.dat

      e <- y.dat - x.dat %*% beta.start
      startSigma <- (t(e) %*% e) / (N.input * min(Time.input) - (K+R))  # OLS for both sigmas

      myPar     <- c(sigma_u = startSigma,
                     sigma_v = startSigma,
                     beta = beta.start,
                     delta = delta.start)

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
        if (length(myPar) != (2 + K + R)){  # 2 sigmas + k betas & r deltas
          stop ("Could not perform estimation.
                A starting point for the estimation must be provided for every parameter.")
        }
        if (!is.numeric(myPar) || myPar[1] <= 0 || myPar[2] <= 0){
          stop ("Could not perform estimation.
                The starting points must be numeric &| sigmas must be <= 0")
        }

        optim.SFM <- nlminb(objective = SFM.within,
                            start = c(myPar),
                            lower = l.int,
                            mu = mu,
                            Time = Time.input,
                            N = N.input,
                            xv = x.dat, y = y.dat, z = z.dat,
                            optim = T)
        }

    } else { # we dont estimate (estimate = F)
      optim.SFM <- NULL
      optim.SFM$par <- myPar  # provided parameters are used as parameters
  }

  # derive the Hessian Matrix based on estimates from the optimization
  hes <- numDeriv::hessian(SFM.within,
                           method = "Richardson",
                           x = optim.SFM$par,
                           xv = x.dat, y = y.dat, z = z.dat,
                           N = N.input,
                           Time = Time.input,
                           mu = mu,
                           optim = T)  # Note optim = T

  # Get & apply variables of the SFM.within model based on estimation  ---------------------------
  # Get all variables from the SFM.within model for further calculations

  if(any(method == c("within", "firstdiff"))){
    if(method == "within"){
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
  }

  if (estimate == F){  # log.ll is obtained from the within model as nlminb is not used.
    optim.SFM$objective <- sum(ret.list$log.ll)
  }

  # assign the estimations
  estimate.sigma_u  <- optim.SFM$par[1]
  estimate.sigma_v  <- optim.SFM$par[2]
  estimate.beta     <- optim.SFM$par[3:(3+K-1)]
  estimate.delta    <- optim.SFM$par[(3+K):(3+R-1)]

  # Calculate Confidence Intervals for the estimates (returns a data frame)
  # alpha can be a vector
  # sigmaCI = NULL & !is.null(sigmaCI) & !is.nan(sigmaCI)
  if ( (!any (sigmaCI <= 0 | sigmaCI > 1)) && !is.null(sigmaCI) && !is.nan(sigmaCI) ){
    conf.Interval <- SFM.CI(estimates = optim.SFM$par, hessianMatrix = hes, alpha = sigmaCI)
  } else if (is.null(sigmaCI)){

  } else if (any (sigmaCI <= 0 | sigmaCI > 1) || is.nan(sigmaCI)){
    cat ("Could not compute Confidence Intervals due to invalid input (sigmaCI must be between [0, 1]")
  }

  # calculate the inefficency index for each panel
  # TODO(Oli): extend inefficency to Time as an vector (check sfm_within for that)
  # inefficency <- SFM.inindex(h = ret.list$h,  # Note h is not within transformed
  #                            sigma2star = ret.list$sigma_2star,
  #                            mu2star = ret.list$mu_2star,
  #                            N = N.input,
  #                            Time = Time.input)

  # Recover the intercept alpha for each panel
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

  # Calculate Model Selection Criterion  ---------------------------
  AIC <- -2 * optim.SFM$objective + 2*length(optim.SFM$par)
  BIC <- -2 * optim.SFM$objective + length(optim.SFM$par) * dim(y.dat)[1]

  # TODO(Clemens): Check AIC / BIC signs

  # Output list
  # TODO(Oli): add option if some values are not created like conf.Interval
  # -> we have currently one fail in unit tester due   to that
  res <- list(call = call, par = myPar, hessian = hes,
              estimates = optim.SFM$par , AIC = AIC, BIC = BIC, estimate = estimate,
              # ci = conf.Interval
              ret.list = ret.list, contrasts = c(attr(myPar, "names")),
              likeihood= optim.SFM$objective)

  class(res) <- c(res$class, "sfmfep")
  res
  # return(list(AIC, BIC))
}

