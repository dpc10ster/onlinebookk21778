rm( list = ls())#mainBinVariability.R
library(caTools);library(ggplot2);source("rocY.R")

mu <- 2;sigma <- 1.5 # experiment with other values
a <- mu/sigma; b <- 1/sigma # a and b parameters
cat("true AUC = ", pnorm(mu/sqrt(1+sigma^2)), "\n")

x <- seq(0.0, 1, 0.01)

zeta <- c(3, 2.5, 2)
FPF <- pnorm(-zeta);FPF <- c(0,FPF,1)
TPF <- pnorm((mu-zeta)/sigma);TPF <- c(0,TPF,1)
pointsData <- data.frame(FPF = FPF, TPF = TPF)
AUC <- trapz(FPF,TPF)
cat("empirical AUC, sparse points = ", AUC, "\n")
rocPlot1 <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = pointsData) + 
  geom_point(data = pointsData)
  
print(rocPlot1)

zeta <- seq(3, -2, -0.5)
FPF <- pnorm(-zeta);FPF <- c(0,FPF,1)
TPF <- pnorm((mu-zeta)/sigma);TPF <- c(0,TPF,1)
pointsData <- data.frame(FPF = FPF, TPF = TPF)
AUC <- trapz(FPF,TPF)
cat("empirical AUC, dense point = ", AUC, "\n")
rocPlot2 <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = pointsData) + 
  geom_point(data = pointsData)
print(rocPlot2)