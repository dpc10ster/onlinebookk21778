#mainBinomialExample1.R
rm(list = ls())
library(exactci)

K2 <- 100;L <- 1;nuP1 <- 0.5;nuP2 <- 0.9;
cat ("K2 = ", K2,
     "\nnuP 1st reader = ", 0.5,
     "\nnuP 2nd reader = ", 0.9,"\n")

seed <- 1;set.seed(seed)
samples1 <- rbinom(K2,L,nuP1)
cat("mean, reader 1 = ", mean(samples1)/L, "\n")

seed <- 1;set.seed(seed)
samples2 <- rbinom(K2,L,nuP2)
cat("mean, reader 2 = ", mean(samples2)/L, "\n")

ret1 <- binom.exact(sum(samples1),K2*L)
ret2 <- binom.exact(sum(samples2),K2*L)

cat ("Rdr. 1: 95% CI = ", ret1$conf.int[1:2],"\n")
cat ("Rdr. 2: 95% CI = ", ret2$conf.int[1:2],"\n")

