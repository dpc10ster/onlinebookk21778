# main4ParameterVs2Parameter.R
rm(list = ls())
source("OpPtsFromZsamples.R")
seed <- 1;set.seed(seed)

# define 4-parameter binormal model
Mu1 <- 3.1;Sigma1 <- 2.3
Mu2 <- 4.7;Sigma2 <- 4
K1 <- 10;K2 <- 12

# get K1/K2 samples from the normal distributions
z1 <- rnorm(K1,mean = Mu1, sd = Sigma1)
z2 <- rnorm(K2,mean = Mu2, sd = Sigma2)

OP <- OpPtsFromZsamples ( z1, z2 )
FPF <- OP$FPF
TPF <- OP$TPF

plot(FPF, TPF)

# if you do not reset the seed, the data will be different
seed <- 1;set.seed(seed)

# define equivalent 2-parameter binormal model
mu <- (Mu2-Mu1)/Sigma1 
sigma <- Sigma2 / Sigma1

# get K1/K2 samples from the normal distributions
z11 <- rnorm(K1) # mu = 0 and sd = 1 is implicit
z22 <- rnorm(K2,mean = mu, sd = sigma)

OP <- OpPtsFromZsamples ( z11, z22 )
FPF1 <- OP$FPF
TPF1 <- OP$TPF

plot(FPF1, TPF1)
