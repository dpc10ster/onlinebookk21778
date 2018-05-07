# MainFrocRaw.R
rm(list = ls())
library(RJafroc)

seed <- 1;set.seed(seed)
mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- -Inf;K1 <- 4;K2 <- 3;
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1))

frocDataRaw <- FROCSimulator(
  mu = mu, lambda = lambda, nu = nu, 
  K1 = K1, K2 = K2, Lk2 = Lk2, zeta1 = zeta1
)

plotFROC <- EmpiricalOpCharac(dataset = frocDataRaw, trts= 1, rdrs = 1, opChType = "FROC")
print(plotFROC$FROCPlot)
