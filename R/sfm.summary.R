summary.sfmfep <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  head <- c("Estimates", "Std.Error", "t-value","confInt")
  mat <- matrix(NA, length(x$contrasts), length(head),
                dimnames = list(x$contrasts, head))
  mat[,1] <- x$estimates; mat[,2] <- x$standerror

  cat("Call:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Estimates:\n")
  printCoefmat(mat)


  cat("\n\nAIC:", paste(deparse(x$aic)))

}


#######-------------tester---------------##########
contrasts <- c("sigma_u","sigma_v","beta","delta")
call <- "y~x1+x2"

cat("\nEstimates:\n")
printCoefmat(mat)


hello <- list(contrasts = c("sigma_u","sigma_v","beta","delta"),
              call = "y~x1+x2",
              standerror = c(1,2,3,4),
              aic= 150,
              estimates=c(1,1,1,1))

summary.sfmfep(x=hello)

