# mainProprocRocs.R
rm( list = ls())
require(caTools)
require(mvtnorm)
require(ggplot2)
source("rocY.R")
source("ProprocFunctions.R")

c1Arr <-   c(-0.1322804, 0.2225588); daArr  <-  c(1.197239,1.740157)
for (i in 1:2)
{
  c1 <- c1Arr[i]
  da <- daArr[i]
  ret <- Transform2ab(da, c1)
  a <- ret$a;b <- ret$b
  if (i == 1) z <- seq(-3, 0, by = 0.01) # may need to adjust limits to view detail of slope plot
  if (i == 2) z <- seq(-3, 5, by = 0.01) # may need to adjust limits to view detail of slope plot
  
  FPF <- seq(0.0, 1, 0.001)
  TPF <- rocY(FPF, a, b)
  
  rocPlot <- data.frame(FPF = FPF, TPF = TPF)
  plotRoc <- ggplot(rocPlot, aes(x = FPF, y = TPF)) + 
    geom_line(size = 2)  + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"))  +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  print(plotRoc)
  
  # do integral numerically
  numAuc <- trapz(FPF, TPF)
  
  # Implement Eqn. 36 from Metz-Pan paper 
  rho <- -(1-c1^2)/(1+c1^2);sigma <- rbind(c(1, rho), c(rho, 1))
  lower <- rep(-Inf,2);upper <- c(-da/sqrt(2),0)
  A_prop <- pnorm(da/sqrt(2)) + 2 * pmvnorm(lower, upper, sigma = sigma)
  A_prop <-  as.numeric(A_prop)
  
  slope <-b*dnorm(a-b*z)/dnorm(-z) # same as likelihood ratio
  
  slopePlot <- data.frame(z = z, slope = slope)
  plotSlope <- ggplot(slopePlot, aes(x = z, y = slope)) + 
    geom_line(size = 2)  + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"))  +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  print(plotSlope)
  
  
  cat("\nC = ", c1, "\nda = ", da, "\nNumerical AUC = ", numAuc, "\nEqn. 36 = ", A_prop,"\n")
}