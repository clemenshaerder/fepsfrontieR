#' @title Plot of the Inefficencys
#' @description Plots the Inefficencys of each panel.
#' @param x is a data frame of the inefficencys per panel.
#' @param ... Additional arguments to the function
#' @importFrom graphics par hist plot abline

plot.sfmfep <- function(x, ...){
  par(mfrow=c(1,2))

  hist(x$Ineff, main = "Hist of Inefficiency", xlab = "Inefficiency")

  plot(y = x$Ineff, x = 1:length(x$Ineff), main = "mean Inefficiency",
       xlab = "panels", ylab = "Inefficiency")
  abline(h = mean(x$Ineff), col= "red")
}
