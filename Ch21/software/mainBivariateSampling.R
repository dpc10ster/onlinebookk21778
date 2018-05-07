# mainBivariateSampling.R
rm(list = ls()) 
library(plotly)
library(MASS);library(mvtnorm)
seed <- 1;set.seed(seed);N <- 10000 #number of samples

mu12  <- 1.5 # mean of diseased distribution, modality 1, truth 2
mu22 <- 2.0  # mean of diseased distribution, modality 2, truth 2
sig1 <- 1.0 # sq. root of variance of diseased distribution, modality 1
sig2 <- 1.5 # sq. root of variance of diseased distribution, modality 2
rho1 <- 0.3 # correlation non-diseased cases
rho2 <- 0.6 # correlation diseased cases
mu1 <- as.vector(array(0, dim = 2)) # vector, non diseased cases
mu2 <- as.vector(c(mu12, mu22)) # vector diseased cases

COV1 <- array(dim = c(2,2))
COV1[1,1] <- 1;COV1[2,2] <- 1;COV1[1,2] <- rho1;COV1[2,1] <- rho1
COV1 <- matrix(COV1, nrow = 2)
z1  <- mvrnorm(N, mu1, COV1)
mu1Est <- c(mean(z1[,1]), mean(z1[,2]))
rho1Est <- cor(z1[,1],z1[,2])
COV1Est <- cov(z1)
cat("expected means of non-diseased cases = ", 0, 0, "\n")
cat("expected correlation of non-diseased cases = ", rho1, "\n")
cat("expected covariance of non-diseased cases = \n")
print(COV1)
cat("\n")
cat("observed means of non-diseased cases = ", mu1Est, "\n")
cat("observed correlation of non-diseased cases = ", rho1Est, "\n")
cat("observed covariance of non-diseased cases = \n")
print(COV1Est)
cat("\n")

COV2 <- array(dim = c(2,2))
COV2[1,1] <- sig1^2;COV2[2,2] <- sig2^2;COV2[1,2] <- 
  rho2*sig1*sig2;COV2[2,1] <- rho2*sig1*sig2
COV2 <- matrix(COV2, nrow = 2)
z2  <- mvrnorm(N, mu2, COV2)
mu2Est <- c(mean(z2[,1]), mean(z2[,2]))
rho2Est <- cor(z2[,1],z2[,2])
COV2Est <- cov(z2)
cat("expected means of diseased cases = ", mu12, mu22, "\n")
cat("expected correlation of diseased cases = ", rho2, "\n")
cat("expected covariance of diseased cases = \n")
print(COV2)
cat("\n")
cat("observed means of diseased cases = ", mu2Est, "\n")
cat("observed correlation of diseased cases = ", rho2Est, "\n")
cat("observed covariance of diseased cases = \n")
print(COV2Est)

upperX <- ceiling(mu12 + 3 * sig1)
upperY <- ceiling(mu22 + 3 * sig2)
x <- seq(-3, upperX, by = 0.1)
y <- seq(-3, upperX, by = 0.1)
z <- array(dim = c(length(x), length(y)))
for (ix in 1:length(x)){
  for (iy in 1:length(y)){
    z[ix, iy] <- max(dmvnorm(c(x[ix], y[iy]), sigma = COV1), 
                     dmvnorm(c(x[ix], y[iy]), mean = c(mu12, mu22), sigma = COV2))
  }
}

p <- plot_ly(x = x, y = y, z = z, type = "surface")
print(p)
