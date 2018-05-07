#mainBinomialExample2.R
rm(list = ls())
library(exactci)

K2 <- c(97,2,1);Lk <- c(1,2,3);nuP1 <- 0.5;nuP2 <- 0.9;
samples1 <- array(dim = c(sum(K2),length(K2)))
cat("K2[1] =", K2[1],
    "\nK2[2] =", K2[2],
    "\nK2[3] =", K2[3], 
    "\nnuP1 =", nuP1, "\nnuP2 =", nuP2, "\n")

seed <- 1;set.seed(seed)
for (l in 1:length(K2)) {
  samples1[1:K2[l],l] <- rbinom(K2[l],Lk[l],nuP1)
}
cat("obsvd. mean, reader 1 = ", 
    sum(samples1[!is.na(samples1)])/sum(K2*Lk), "\n")

samples2 <- array(dim = c(sum(K2),length(K2)))
seed <- 1;set.seed(seed)
for (l in 1:length(K2)) {
  samples2[1:K2[l],l] <- rbinom(K2[l],Lk[l],nuP2)
}
cat("obsvd. mean, reader 2 = ", 
    sum(samples2[!is.na(samples2)])/sum(K2*Lk), "\n")

ret1 <- binom.exact(sum(samples1[!is.na(samples1)]),sum(K2*Lk))
ret2 <- binom.exact(sum(samples2[!is.na(samples2)]),sum(K2*Lk))

cat ("Rdr. 1: 95% CI = ", ret1$conf.int[1:2],"\n")
cat ("Rdr. 2: 95% CI = ", ret2$conf.int[1:2],"\n")