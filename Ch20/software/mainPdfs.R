# mainPdfs.R
rm( list = ls())
require(ggplot2)
source("rocY.R")

aArray <- c(0.7, 0.7, 1.5, 2)
bArray <- c(0.5, 1.5, 0.5, 0.5)
z1 <- seq(-5, 3, by = 0.01)
z2 <- seq(-5, 7, by = 0.01)
FPF <- seq(0.0, 1, 0.01)

pdf1 <- dnorm(z1)

for (i in 1:length(aArray))
{
  a <- aArray[i]
  b <- bArray[i]
  TPF <- rocY(FPF, a, b)
  rocPlot <- data.frame(FPF = FPF, TPF = TPF)
  plotRoc <- ggplot(
    rocPlot, aes(x = FPF, y = TPF)) + 
    geom_line()
  print(plotRoc)
  
  pdf2 <- dnorm(z2, a/b, sd = 1/b)
  df <- data.frame(
    z = c(z1, z2), pdfs = c(pdf1, pdf2), 
    truth = c(rep('non-diseased', length(pdf1)), 
              rep('diseased', length(pdf2))))
  
  rocpdfs <- ggplot(
    df, 
    aes(x = z, y = pdfs, color = truth)) + 
    geom_line()
  print(rocpdfs)
  
  cat("a = ", a, ", b = ", b, "\n")
  break
}