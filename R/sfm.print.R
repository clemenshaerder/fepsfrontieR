#' @title Print Output for sfmfep
#' @description A print of valuable variables is done when sfmfep is called.
#' @param ... Additional arguments to the function
#' @param x is a data frame of the inefficencys per panel.
#' @param digits defines the amount of digits for the values of the summary
#' @export

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

