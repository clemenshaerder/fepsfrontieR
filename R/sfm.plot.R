plot.sfmfep <- function(x, ...){
  par(mfrow=c(1,2))

  hist(x$Ineff, main = "Hist of Inefficiency", xlab = "Inefficiency")

  plot(x$Ineff, main = "mean Inefficiency",
       xlab = "panels", ylab = "Inefficiency")
  abline(h = mean(x$Ineff), col= "red")
}

