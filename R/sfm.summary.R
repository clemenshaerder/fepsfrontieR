summary.sfmfep <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if(x$estimate){

    head1 <- c("Coefficients", "Std.Error", "T-stat")
    head2 <- c("alpha", "inefficiency")

    mat1 <- matrix(NA, length(x$contrasts), length(head1),
                   dimnames = list(x$contrasts, head1))
    mat2 <- matrix(NA, length(x$alpha), length(head2),
                   dimnames = list(rownames(x$alpha), head2))

    T.stat <- x$coefficients / x$standerror
    T.stat[c(1,2)] <- NA

    mat1[, 1] <- x$coefficients
    mat1[, 2] <- x$standerror
    mat1[, 3] <- T.stat
    mat1 <- cbind(mat1, x$conf)

    mat2[, 1] <- x$alpha
    mat2[, 2] <- x$Ineff

    cat("Call:\n",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    if(x$bootstrap == T){
      cat("method:", deparse(x$method), "\nEstimates are based on ", x$B , "bootstrap samples\n\n")
    } else {
      cat("method:", deparse(x$method), "\n\n")
    }

    if(!is.null(x$conf)){
      cat("Coefficients:\t\t\t\t\tconfInterval\n")
    } else {
      cat("Coefficients:\n")
    }

    printCoefmat(mat1)

    cat("\nlog.Likelihood:", deparse(round(x$objective, digits)),
        "; AIC:", deparse(round(x$aic, digits)),
        "; BIC:", deparse(round(x$bic, digits)), "\n\n")

    printCoefmat(mat2)
  } else {
    cat("Likelihood:", deparse(round(x$objective, digits)))
 }
}

