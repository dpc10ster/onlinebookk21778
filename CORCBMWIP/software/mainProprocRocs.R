# mainProprocRocs.R
rm( list = ls())
require(ggplot2)

cArr <-   c(0.9999, 0.9999); daArr  <-  c(0, 0)
b <- -(cArr + 1)/(cArr - 1)
a <- (daArr/sqrt(2)) * sqrt(1 + b^2)
for (i in 1:2)
{
  z <- seq(-3, 5, by = 0.01) 
  
  FPF <- seq(0.0, 1, 0.001)
  TPF <- pnorm(a[i] + b[i] * qnorm(FPF))
  
  rocPlot <- data.frame(FPF = FPF, TPF = TPF)
  plotRoc <- ggplot(rocPlot, aes(x = FPF, y = TPF)) + geom_line()
  print(plotRoc)
  
}
