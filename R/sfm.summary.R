#' @title Summarizes a sfmfep fitted model
#' @description  summary.sfmfep is the sfmfep specific method for
#' the generic function summary which summarize objects
#' returned by modelling functions.
#' @param ... Additional arguments to the function
#' @param object is a data frame of the inefficencys per panel.
#' @param digits defines the amount of digits for the values of the summary
#' @export

summary.sfmfep <- function(object, digits = max(3L, getOption("digits") - 3L), ...){

  if(object$estimate){

    head2 <- c("alpha", "inefficiency")

    mat2 <- matrix(NA, length(object$alpha), length(head2),
                   dimnames = list(rownames(object$alpha), head2))

    mat2[, 1] <- object$alpha
    mat2[, 2] <- object$Ineff

    if(!is.null(object$standerror)){

      head1 <- c("Coefficients", "Std.Error")

      mat1 <- matrix(NA, length(object$contrasts), length(head1),
                   dimnames = list(object$contrasts, head1))

       mat1[, 1] <- object$coefficients
       mat1[, 2] <- object$standerror

        if(!is.null(object$conf)){
          mat1 <- cbind(mat1, object$conf)
        }
       } else {

         head1 <- c("Coefficients")

         mat1 <- matrix(NA, length(object$contrasts), length(head1),
                        dimnames = list(object$contrasts, head1))

         mat1[, 1] <- object$coefficients
         if(!is.null(object$conf)){
           mat1 <- cbind(mat1, object$conf)
         }
      }

    cat("Call:\n",
        paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

    if(object$bootstrap == T){
      cat("method:", deparse(object$method), "\nEstimates are based on ",
          object$B , "bootstrap samples\n\n")
    } else {
      cat("method:", deparse(object$method), "\n\n")
    }

    if(!is.null(object$conf)){
      cat("Coefficients:\t\t\t\t\tconfInterval\n")
    } else {
      cat("Coefficients:\n")
    }

    printCoefmat (mat1)

    cat("\nlog.Likelihood:", deparse(round(object$objective, digits)),
        "; AIC:", deparse(round(object$aic, digits)),
        "; BIC:", deparse(round(object$bic, digits)), "\n\n")

    printCoefmat (mat2)
  } else {
    cat("Likelihood:", deparse(round(object$objective, digits)))
 }
}

