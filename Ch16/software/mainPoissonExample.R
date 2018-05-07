#mainPoissonExample.R
rm(list = ls())
library(exactci)

K <- 100
lambdaP1 <- 1;lambdaP2 <- 2;
cat ("K = ", K, 
     "\nlambdaP 1st reader = ", lambdaP1,
     "\nlambdaP 2nd reader = ", lambdaP2,"\n")

seed <- 1;set.seed(seed)
samples1 <- rpois(K,lambda = lambdaP1)
cat("obs. mean, reader 1 = ", mean(samples1), "\n")

seed <- 1;set.seed(seed)
samples2 <- rpois(K,lambda = lambdaP2)
cat("obs. mean, reader 2 = ", mean(samples2), "\n")

ret1 <- poisson.exact(sum(samples1),K)
ret2 <- poisson.exact(sum(samples2),K)

cat ("Rdr. 1: 95% CI = ", ret1$conf.int[1:2],"\n")
cat ("Rdr. 2: 95% CI = ", ret2$conf.int[1:2],"\n")
