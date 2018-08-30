print.sfmfep <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if(x$estimate){
    x$out <- setNames(x$coefficients, x$contrasts)

    cat("\nCall: ", "\n",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Coefficients")

    cat(":\n")
    print.default(format(x$out, digits = digits),
                  print.gap = 2, quote = FALSE)
  } else {
    cat("Likelihood:", deparse(round(x$objective, digits)))
  }

}

