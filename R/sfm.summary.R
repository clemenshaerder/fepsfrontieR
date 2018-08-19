summary.sfmfep <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  head1 <- c("Estimates", "Std.Error", "t-value")
  head2 <- c("alpha", "inefficiency")

  mat1 <- matrix(NA, length(x$contrasts), length(head1),
                dimnames = list(x$contrasts, head1))
  mat2 <- matrix(NA, length(x$alpha), length(head2),
                 dimnames = list(1:length(x$alpha), head2))


  mat1[, 1] <- x$estimates
  mat1[, 2] <- x$standerror
  mat1[, 3] <- x$tvalue
  mat1 <- cbind(mat1, x$conf)

  mat2[, 1] <- x$alpha
  mat2[, 2] <- x$Ineff

  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if(x$bootstrap == T){
    cat("Estimates are based on ", x$B , "bootstrap samples\n")
  } else {
    cat("method:", deparse(x$method), "\n")
  }

  cat("Estimates:\t\t\t\tconfInterval\n")
  printCoefmat(mat1)

  cat("log.Likelihood:", deparse(round(x$objective, digits)),
      "; AIC:", deparse(x$aic),
      "; BIC:", deparse(x$bic), "\n\n")

  printCoefmat(mat2)

}

#######-------------tester---------------##########
#contrasts <- c("sigma_u","sigma_v","beta","delta")
#call <- "y~x1+x2"
#cat("\nEstimates:\n")
#printCoefmat(mat)


#hello <- list(contrasts = c("sigma_u","sigma_v","beta1" ,"beta2", "delta1", "delta2"),
#              call = "y~x1+x2",
#              standerror = c(1,2,3,4,5,6),
#              aic = 150,
#              estimates = c(1,1,1,1,1,1),
#              bootstrap = T,
#              tvalue = c(1,1,1,1,1,1),
#              conf = h1$conf.Interval,
#              objective = h1$objective,
#              bic = 140,
#              alpha = 1:10,
#              Ineff = rnorm(10, 0, 1),
#              method = "first Diff",
#              B = 1000
#              )

#summary.sfmfep(x=hello)




