rm(list = ls()) # MainBivariateDensity.R
library(mvtnorm)

X <- c(0.1,0.2)
mu1 <- 0;mu2 <- 0
rho <- 0.5;var1 <- 1.5;var2 <- 1.5;
COV <- matrix(c(var1, rho*sqrt(var1*var2), rho*sqrt(var1*var2), var2),2)
ret  <- dmvnorm(X, sigma = COV)
cat("rho =", rho, ", var1 =", var1, ", var2 =", var2, "\n")
cat("density at 0.1, 0.2 = ", ret, "\n")