#' @title Fixed-Effect Stochastic Frontier Model Estimation for Panel Data
#' @description Estimates a Stochastic Frontier Model for Fixed-Effects.
#'     sfmfep is used to fit fixed-effect stochastic frontier models for panel data
#'     using a specified model transformation. Bootstrapping can be performed to
#'     calculate the standard errors instead of a numerical deriviation
#'     via the hessian matrix.
#' @param formula an object of class "formula" in the form of
#'     y ~ x1 + ... + x_k + (z1 + ... + z_r)
#' @param data an optional data frame, list or environment (or object coercible by
#'     as.data.frame to a data frame) containing the variables in the model.
#' @param panel an optional vector specifying the panels to be used in the fitting process.
#' @param N an optional integer specifying the total amount of panels in the data set.
#'     If N is entered, Time is also a required input.
#' @param Time an optional integer specifying the amount of observations per panel.
#'     If Time is entered, N is also a required input.
#' @param method a required string specifying the method ("within" or "firstdiff").
#' @param mu is the mean of a truncated normal distribution of the stochastic inefficiency.
#' @param sigmaCI is an optional vector specifying the significance
#'     values of the confidence intervals for the MLE estimates.
#' @param estimate TRUE or FALSE specifies if "myPar" is used as
#'     starting point of the estimation, or if "myPar" is used to fit
#'     a given stochastic frontier model.
#' @param bootstrap is an optional boolean variable. If it is set to TRUE, bootstrapping is
#'     performed. "B" must be specified.
#' @param B is an integer which is a required input for Bootstrap = T. It defines
#'     the amount of bootstrap samples.
#' @param parallel is an optional boolean variable. If it is set to TRUE, bootstrapping is
#'     performed with parallelization, using all available cores - 1.
#'     Only available for OS Windows.
#' @param myPar is a vecor which has to be entered in the following order:
#'     c(sigma_v, sigma_u, beta = c( ), delta = c( ))
#' @return
#'     sfmfep returns an object of class S3 "sfmfep". The function summary( )
#'     can be used to obtain or print a summary of the results.
#' @examples
#' # Fit of a simple model with balanced panels.
#' # Definition with *N & T*.
#' # Data has 10 observations for each of the 20 panels.
#'
#' fit1 <- sfmfep(formula = y ~ x + (z),
#'     method = "within", N = 20, Time = 10, data = sfm.data)
#' summary(fit1)
#'
#'# ---------------
#' # Fit of a simple model with balanced panels.
#' # Definition with *panels*.
#'
#' fit2 <- sfmfep(formula = y ~ x + (z),
#'      method = "within", panel = sfm.data$producer, data = sfm.data)
#' summary(fit2)
#'
#' # ---------------
#' # Fit of a simple model with *Bootstrapping* using
#' # *method = firstdiff*
#' # with different *sigmas* for *Confidence Intervals*
#'
#' fit3 <- sfmfep(formula = y ~ x + (z),
#'     bootstrap = TRUE, B = 20, sigmaCI = c(0.1, 0.05),
#'     method = "firstdiff", panel = sfm.data$producer, data = sfm.data)
#' summary(fit3)
#'
#' # ---------------
#' # Fitting a model *without estimating*, providing model parameters.
#'
#' fit4 <- sfmfep(formula = y ~ x + (z),
#'     myPar = c(sigma_u = 0.2, sigma_v = 0.1, beta = 0.5, delta = 0.5),
#'     estimate = FALSE,
#'     method = "firstdiff", panel = sfm.data$producer, data = sfm.data)
#' summary(fit4)
#'
#' # ---------------
#' # Perform Bootstrapping with starting points for the optimizier
#'
#' fit5 <- sfmfep(formula = y ~ x + (z),
#'     myPar = c(sigma_u = 0.2, sigma_v = 0.1, beta = 0.5, delta = 0.5),
#'     estimate = TRUE, bootstrap = TRUE, B = 20,
#'     method = "firstdiff", panel = sfm.data$producer, data = sfm.data)
#' summary(fit5)
#'
#' @importFrom magrittr %>%
#' @importFrom numDeriv hessian
#' @importFrom stringr str_extract
#' @importFrom numDeriv hessian
#' @importFrom dplyr arrange_
#' @importFrom stats as.formula dnorm nlminb pnorm rnorm runif sd setNames printCoefmat
#' @importFrom MASS ginv
#' @export

sfmfep <- function(formula, data, panel = NULL, N = NULL, Time = NULL,
                   method = "firstdiff", mu = 0,  sigmaCI = 0.05, estimate = T,
                   bootstrap = F, B = NULL, parallel = F, myPar = c(sigma_u = NULL,
                   sigma_v = NULL, beta = c(NULL), delta = c(NULL))){

  call <- match.call ()

  # Error handling of input ---------------------------

  # Tests if a correct formula has been plugged in (i.e. y ~ x...)
  # We do not include the option to add y ~ . , as it is not possible to identify the z variables.
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

  # Check the format of the data input
  if ((is.data.frame (data) || is.matrix (data)) == F){
    stop ("Your data should be a data.frame or matrix")
  }

  # Tests if any of the model variables are not numeric
  if (any (sapply (sel.data, function (sel.data) !is.numeric (sel.data)))){
    stop ("All variables must be numeric.")
  }

  # Tests if any NAs exist in the data sel.data
  dataInvalid <- sapply (sel.data, function (sel.data) sum (is.na (sel.data)) +
                                                       sum (is.infinite (sel.data)) +
                                                       sum (is.null (sel.data)))
  if (sum (dataInvalid) > 0){
    print ("Function Stopped. You have NAs/Inf/NULL in your selected data. Summary:")
    stop (print (dataInvalid))
  }

  # Test if method is correctly specified
  if (!any (method == c("within", "firstdiff"))){
    stop ("Invalid input. A valid *method* must be chosen *within* or *firstdiff*")
  }

  # Test if mu is correctly specified
  if (mu < 0 || is.na(mu) == T || is.null(mu) == T){
    stop ("The mean of the truncated normal distribution (*mu*) must be a positive value")
  }

  # Tests if the sig.niveau sigmaCI for the Confidence Interval are correctly entered
  if (!is.null (sigmaCI)){ # NULL would throw an error, but is a valid input.
    if (any (sigmaCI <= 0 | sigmaCI > 1) || is.nan (sigmaCI)){
      cat ("Can not compute Confidence Intervals due to invalid input (sigmaCI must be between [0, 1]")
    }
  }

  # Test if estimate is correctly specified
  if (!any (estimate == c(T, F))){ # T == TRUE
    stop ("Invalid input. *estimate* must be either True or False")
  }

  # Test if bootstrap is correctly specified
  if (!any (bootstrap == c(T, F))){
    stop ("Invalid input. *bootstrap* must be either True or False.")
  }

  # Tests if B is correctly specified if bootstrapping is performed
  if (bootstrap == T){
    if (B <= 0 || !is.numeric(B) || length(B) > 1){
      stop("Could not perform calculations.
           B must be an integer > 0")
    }
  }

  # Tests if parallel is correctly specified if bootstrapping is performed
  if (!any (parallel == c(T, F))){
    stop ("Invalid input. *parallel* must be either True or False")
  }


  # Data Wrangling & Error handling of group, N & T  & myPar ---------------------------

  formula <- as.character (formula)

  # extracts z-variables (in curly brackets) by regular expressions
  zForm <- str_extract (formula, "(?<=\\().*?(?=\\))")[3]
  # create a new formula only with z variables to continue data wrangling
  extractZ <- as.formula (paste (formula[2], formula[1], zForm))

  R <- length (all.vars (extractZ)) - 1  # -1, as extractZ is of form *y* ~ z1 + ... + zr
  if (R <= 0){
    stop ("No required z variable were defined in the formula.
          Please add at least one z variable in your formula in brackets (z_i)")
  }

  totCountVar <- dim (sel.data)[2]  # total amount of variables
  K <- totCountVar - R - 1 # all variables - r Z-variables - (1) y-variable = K-variables

  # Tests if myPar is correctly specified if it is not NULL
  if (!is.null(myPar)){
    if (length (myPar) != (2 + K + R)){  # 2 sigmas + k betas & r deltas
      stop ("Could not perform estimation.
            A starting point for the estimation must be provided for every parameter.")
    }
    if (!is.numeric (myPar) || myPar[1] <= 0 || myPar[2] <= 0){
      stop ("Could not perform estimation.
            The starting points must be numeric and/or sigmas must be <= 0")
    }
  }

  # N & T, or panel must be assign, else we can not compute N & T
  # If only panel and N & T is assigned, we use panel
  # First, we check if no option is defined
  if (length (panel) > 1 & is.list (panel) == F){  # If panel is defined as panel$var
    # get the column of the specified vector (e.g. data$panel)
    try (colPanel <- which(data == panel,  arr.ind = T)[1,2], silent = F)
    if (!exists("colPanel")){
      stop ("*panel* must be a column in your data.")
    }
    panel <- colnames (data)[colPanel]
  } else if (is.list (panel)){  # if panel is defined as panel[, col]
    panel <- colnames (panel)
  }
  if (length(panel) > 1){
    stop ("Incorrect defintion of *panel*. Could not perform calculations.")
  }
  if ((is.null (N) | is.null (Time)) & is.null (panel) ){
    stop ("You have to either specify N = panels & Time = obs. per panel
          or provide a panel column")

  # Second, we check if N & T are both well definded and if panel is not
  } else if ((is.double (N) && is.double (Time)) && !is.character (panel)){

    N.input <- N; Time.input <- Time

    # We check if N & T match the data dimensions
    if (length (Time.input) == 1){
      if ( sum (N.input * Time.input) != dim (sel.data)[1]) {
        stop ("Your input (N & T) doesn not match the data dimensions.
              Calculations could not be performed.")
      }
    } else {
      if ( sum (Time.input) != dim (sel.data)[1]) {
        stop ("Your input (N & T) doesn not match the data dimensions.
            Calculations could not be performed.")
      }
    }
  # Third step: if panel are specified, we check if the panels exist as column name
  } else {
    if (is.null (panel) || try (exists (panel, data) == F, silent = T)){
      stop ("Couldnt match input *panel* with colnames")
    } else {
      # Finally, if only panel is chosen or all three options we use panel to get N & T
      sel.data <- cbind (sel.data, data[panel])  # adds panel to the selected data
      sel.data <- sel.data %>% arrange_(.dots = panel)  # sort the panels
      N.input <- dim (table (sel.data[panel]))[1]
      Time.input <- table (sel.data[panel])

      # get names of panels for proper output
      panelName <- unique (sel.data[panel]) %>% arrange_(.dots = panel)
    }
  }

  # Adjustment for balanced panels, which is required for other functions
  if(length (Time.input) == 1){
    Time.input <- rep (Time, N)
  }

  # cumTime servves as an index for further computations in "estimation"
  cumTime <- c(0, cumsum (Time.input))

  y.dat <- as.matrix(sel.data[, 1])  # y is always the first value in the formula
  x.dat <- as.matrix(sel.data[, 2:(1+K)])
  z.dat <- as.matrix(sel.data[, (1+K+1):(1+K+R)])


  # Estimation  ---------------------------

  #  lower boundary for optimization
  l.int <- c(0.0001, 0.0001, rep(-Inf, K), rep(-Inf, R))  # Variation can not be negative

  # Program flow of esimtation ("->" stands for a performed action)
  # estimate T
  #   myPar = NULL -> calculate starting Points for optimzier
  #     bootstrap F -> estimate with within / firstDiff
  #     bootstrap T -> bootstrap with within / firstDiff
  #   myPar = ...
  #     bootstrap F -> estimate with within / firstDiff
  #     bootstrap T -> bootstrap with within / firstDiff
  # estimate F
  #   fit with within or firstDiff

  if (estimate == T){
    if (is.null(myPar) == T){  # we generate appropriate starting values

      beta.start  <- solve (t(x.dat) %*% x.dat) %*% t(x.dat) %*% y.dat  # OLS for beta
      delta.start <- solve (t(z.dat) %*% z.dat) %*% t(z.dat) %*% y.dat  # OLS for delta

      # OLS for both sigmas
      e           <- y.dat - x.dat %*% beta.start
      startSigma  <- (t(e) %*% e) / (sum (Time.input) - (K + R))

      myPar       <- c(sigma_u = startSigma,
                      sigma_v = startSigma,
                      beta = beta.start,
                      delta = delta.start)

      if (bootstrap == F) {
        if (method == "within"){
          optim.SFM <- nlminb (objective = SFM.within,
                               lower = l.int,
                               start = myPar,
                               Time = Time.input,
                               N = N.input,
                               xv = x.dat,
                               y = y.dat,
                               z = z.dat,
                               mu = mu,
                               K = K,
                               R = R,
                               cumTime = cumTime,
                               optim = T)
        } else {
          optim.SFM <- nlminb (objective = SFM.firstDiff,
                               lower = l.int,
                               start = myPar,
                               Time = Time.input,
                               N = N.input,
                               xv = x.dat,
                               y = y.dat,
                               z = z.dat,
                               mu = mu,
                               K = K,
                               R = R,
                               cumTime = cumTime,
                               optim = T)
        }
      } else { # else bootstrap = T use bootstrapping
        optim.SFM <- SFM.bootstrap (y = y.dat,
                                    xv = x.dat,
                                    z = z.dat,
                                    mu = mu,
                                    N = N.input,
                                    Time = Time.input,
                                    myPar = myPar,
                                    R = R,
                                    K = K,
                                    B = B,
                                    method = method,
                                    lowerInt = l.int,
                                    cumTime = cumTime,
                                    sigmaCI = sigmaCI,
                                    parallel = parallel)
      }
    } else { # else myPar is defined & not NULL
      if (bootstrap == F){
        if (method == "within"){
          optim.SFM <- nlminb (objective = SFM.within,
                               start = c(myPar),
                               lower = l.int,
                               mu = mu,
                               Time = Time.input,
                               N = N.input,
                               cumTime = cumTime,
                               K = K,
                               R = R,
                               xv = x.dat,
                               y = y.dat,
                               z = z.dat,
                               optim = T)
        } else {  # else use first-difference method
          optim.SFM <- nlminb (objective = SFM.firstDiff,
                               start = c(myPar),
                               lower = l.int,
                               mu = mu,
                               Time = Time.input,
                               N = N.input,
                               cumTime = cumTime,
                               K = K, R = R,
                               xv = x.dat,
                               y = y.dat,
                               z = z.dat,
                               optim = T)
        }
      } else {  # else bootstrap = T use bootstrapping
        optim.SFM <- SFM.bootstrap (y = y.dat,
                                    xv = x.dat,
                                    z = z.dat,
                                    mu = mu,
                                    N = N.input,
                                    Time = Time.input,
                                    myPar = myPar,
                                    method = method,
                                    R = R, K = K,
                                    B = B,
                                    cumTime = cumTime,
                                    lowerInt = l.int,
                                    sigmaCI = sigmaCI,
                                    parallel = parallel)
      }
    }
  } else { # we dont estimate (estimate = F)
      # As we dont estimate optim.SFM list doesn not exit.
      # To perform further calculations we create and assign essential inputs.
      optim.SFM     <- NULL
      if (is.null (myPar)){
        stop ("*myPar* must be specified to fit your model.")
      }
      optim.SFM$par <- myPar  # provided parameters are used as parameters
  }


  # Fit the model based on the estimation  ---------------------------

  # Fit with the specified model
  if (method == "within"){
    ret.list <- SFM.within (optim = F,  # Note optim = F
                            N = N.input,
                            Time = Time.input,
                            xv = x.dat , y = y.dat, z = z.dat,
                            K = K, R = R,
                            cumTime = cumTime,
                            par = optim.SFM$par)
  } else {  # else use firstDiff
    ret.list <- SFM.firstDiff (optim = F,  # Note optim = F
                               N = N.input,
                               Time = Time.input,
                               xv = x.dat, y = y.dat, z = z.dat,
                               cumTime = cumTime,
                               K = K, R = R,
                               par = optim.SFM$par)
  }

  # log.ll is obtained from the model as nlminb is not primarily used.
  if (estimate == F || bootstrap == T){
    optim.SFM$objective <- sum (ret.list$log.ll * -1)
  }

  # Need an error handler when the model isnt well defined
  if (is.nan (optim.SFM$objective) || optim.SFM$objective == Inf ||
      optim.SFM$objective < 0){
    stop ("Optimizer nlminb( ) could not find a valid solution.
           Try different starting points with *myPar* or adapt the model.")
  }


  # Calculate Inefficencys  ---------------------------

  inefficency <- SFM.inindex (h = ret.list$h, sigma2star = ret.list$sigma_2star,
                              cumTime = cumTime, mu2star = ret.list$mu_2star,
                              N = N.input)
  if (!is.null (panel)){
    # if panel is specified, we assign the respective names to the inefficencys
    rownames (inefficency) <- as.matrix (panelName)
  }


  # Recover Alpha (fixed-effects) ---------------------------

  estimate.sigma_u  <- optim.SFM$par[1]
  estimate.sigma_v  <- optim.SFM$par[2]
  estimate.beta     <- optim.SFM$par[3:(3+K-1)]

  alpha <- SFM.alpha (y = y.dat , x = x.dat, beta = estimate.beta,
                      sigma_u = estimate.sigma_u, sigma_v = estimate.sigma_v,
                      h = ret.list$h.trans, epsilon = ret.list$eps.trans,
                      cumTime = cumTime,
                      N = N.input, Time = Time.input, mu = mu)

  if (!is.null (panel)){
    # if panel is specified, we assign the respective names to the inefficencys
    rownames(alpha) <- as.matrix(panelName)
  }


  # Calculate Confidence Intervals  ---------------------------

  # derive the Hessian Matrix based on estimates from the optimization
  if (bootstrap == F){ # not required when bootstrapping is chosen
    if (method == "within"){
      hes <- hessian(SFM.within,
                     x = optim.SFM$par,
                     xv = x.dat,
                     y = y.dat,
                     z = z.dat,
                     N = N.input,
                     Time = Time.input,
                     mu = mu,
                     cumTime = cumTime,
                     K = K,
                     R = R,
                     # optim = T required for computation
                     optim = T)
    } else { # else use firstDiff
      hes <- hessian(SFM.firstDiff,
                     x = optim.SFM$par,
                     xv = x.dat,
                     y = y.dat,
                     z = z.dat,
                     N = N.input,
                     Time = Time.input,
                     mu = mu,
                     cumTime = cumTime,
                     K = K,
                     R = R,
                     # optim = T required for computation
                     optim = T)
    }
  }

  df = sum (Time.input) - length (optim.SFM$par)  # get degrees of freedom

  # sigmaCI can be a vector of significance levels
  if (bootstrap == F){
    if (!is.null (sigmaCI)){
      # A data frame is returned
      c.Interval <- SFM.CI (estimates = optim.SFM$par, hessianMatrix = hes,
                            alpha = sigmaCI, df = df)
      if(!is.character (c.Interval)){  # split the return for the output
        conf.Interval <- c.Interval[, c(2, 3)]
        standerror    <- c.Interval$standerror
      } else {  # if no CI is returned
        conf.Interval <- NULL
        standerror <- NULL
      }
      # If the optimizer is not finding a valid ouput hessian(..)
      # creates a matrix of nans # check if this is still required
    } else if( any (is.nan (hes) == T)){
      conf.Interval <- as.matrix(rep(NA, length(optim.SFM$par)))
      standerror    <- as.matrix(rep(NA, length(optim.SFM$par)))
    } else { # if no CI is wanted we return NAs
      conf.Interval <- as.matrix(rep(NA, length(optim.SFM$par)))
      standerror    <- as.matrix(rep(NA, length(optim.SFM$par)))
    }
  } else { # else bootstrap = T
    conf.Interval <- optim.SFM$conf.Interval # CI is calculated by Bootstrapping
    standerror    <- optim.SFM$standerror
  }


  # Model Selection Criterion  ---------------------------

  # Adjusted formula with "* -1" as the used objective is the "-"log.ll (needed for optimization)
  AIC <- -2 * -1 * optim.SFM$objective + 2 * length (optim.SFM$par)
  BIC <- -2 * -1 * optim.SFM$objective + length (optim.SFM$par) * log (dim (y.dat)[1])


  # Output  ---------------------------

  if(estimate == F){
    res <- list (objective = optim.SFM$objective,
                 estimate = estimate)
  } else {
    res <- list (call = call,
                 par = myPar,
                 coefficients = optim.SFM$par,
                 aic = AIC,
                 bic = BIC,
                 estimate = estimate,
                 conf = conf.Interval,
                 alpha = alpha,
                 Ineff = inefficency,
                 # hessian = hes,
                 standerror = standerror,
                 contrasts = c(attr (optim.SFM$par, "names")),
                 objective = optim.SFM$objective,
                 bootstrap = bootstrap,
                 method = method,
                 B = B,
                 tvalue = NULL)
  }

  class(res) <- c(res$class, "sfmfep")
  res
}


