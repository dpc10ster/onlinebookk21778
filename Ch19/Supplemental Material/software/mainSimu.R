rm(list = ls()) # mainRSM.R 
library(RJafroc)

seed <- 1
set.seed(seed)
K1 <- 5000;K2 <- 5000
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1))
mu <- 2;lambda <- 1.5;nu <- 1 ;zeta1 <- -2
simuDataRaw <- SimulateFrocData(
  mu = mu, lambda = lambda, nu = nu, 
  K1 = K1, K2 = K2, Lk2 = Lk2, zeta1 = zeta1
)

simuData <- DfBinDataset(simuDataRaw, 6)

retRsm <- FitRsmRoc(simuData, 1, 1)

