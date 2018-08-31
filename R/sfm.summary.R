#' @title Summarizes a sfmfep fitted model
#' @description  summary.sfmfep is the sfmfep specific method for
#' the generic function summary which summarize objects
#' returned by modelling functions.
#' @param ... Additional arguments to the function
#' @param x is a data frame of the inefficencys per panel.
#' @param digits defines the amount of digits for the values of the summary
#' @importFrom stats printCoefmat

summary.sfmfep <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if(x$estimate){

    head2 <- c("alpha", "inefficiency")

    mat2 <- matrix(NA, length(x$alpha), length(head2),
                   dimnames = list(rownames(x$alpha), head2))

    mat2[, 1] <- x$alpha
    mat2[, 2] <- x$Ineff

    if(!is.null(x$standerror)){

      head1 <- c("Coefficients", "Std.Error")

      mat1 <- matrix(NA, length(x$contrasts), length(head1),
                   dimnames = list(x$contrasts, head1))

       mat1[, 1] <- x$coefficients
       mat1[, 2] <- x$standerror

        if(!is.null(x$conf)){
          mat1 <- cbind(mat1, x$conf)
        }
       } else {

         head1 <- c("Coefficients")

         mat1 <- matrix(NA, length(x$contrasts), length(head1),
                        dimnames = list(x$contrasts, head1))

         mat1[, 1] <- x$coefficients
         if(!is.null(x$conf)){
           mat1 <- cbind(mat1, x$conf)
         }
      }

    cat("Call:\n",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    if(x$bootstrap == T){
      cat("method:", deparse(x$method), "\nEstimates are based on ",
          x$B , "bootstrap samples\n\n")
    } else {
      cat("method:", deparse(x$method), "\n\n")
    }

    if(!is.null(x$conf)){
      cat("Coefficients:\t\t\t\t\tconfInterval\n")
    } else {
      cat("Coefficients:\n")
    }

    printCoefmat (mat1)

    cat("\nlog.Likelihood:", deparse(round(x$objective, digits)),
        "; AIC:", deparse(round(x$aic, digits)),
        "; BIC:", deparse(round(x$bic, digits)), "\n\n")

    printCoefmat (mat2)
  } else {
    cat("Likelihood:", deparse(round(x$objective, digits)))
 }
}

