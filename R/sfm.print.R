print.sfmfep <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if(x$estimate == T){
    x$out <- setNames(x$estimates, x$contrasts)

    cat("\nCall: ", "\n",
        paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("Estimates")

    cat(":\n")
    print.default(format(x$out, digits = digits),
                  print.gap = 2, quote = FALSE)
  } else {
    paste(deparse(x$likelihood))
  }

}
