# mainImproperRocs.R
rm( list = ls())
require(ggplot2)
source("rocY.R")

aArray <- c(0.7, 0.7, 1.5, 2)
bArray <- c(0.5, 1.5, 0.5, 0.5)
z <- seq(-3, 5, by = 0.01) 
FPF <- seq(0.0, 1, 0.01)
for (i in 1:length(aArray))
{
  a <- aArray[i]
  b <- bArray[i]
  TPF <- rocY(FPF, a, b)
  rocPlot <- data.frame(FPF = FPF, TPF = TPF)
  p <- ggplot(rocPlot, aes(x = FPF, y = TPF)) + 
    geom_line(size = 2) + 
    scale_color_manual(values = "black") + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"))  +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  print(p)
  cat("a = ", a, ", b = ", b, "\n")
}
