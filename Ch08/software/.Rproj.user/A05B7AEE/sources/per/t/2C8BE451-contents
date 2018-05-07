rm(list = ls()) # mainTypeIErrors.R
source("Wilcoxon.R")

seed <- 1;set.seed(seed)
mu <- 1.5;sigma <- 1.3;K1 <- 50;K2 <- 52

# cheat to find the population mean and std. dev.
AUC <- array(dim = 10000)
for (i in 1:length(AUC)) {
  zk1 <- rnorm(K1);zk2 <- rnorm(K2, mean = mu, sd = sigma)  
  AUC[i] <- Wilcoxon(zk1, zk2)
}
sigmaAUC  <-  sqrt(var(AUC));muAUC   <-  mean(AUC)

nTrials <- 2000
alpha <- 0.05 # size of test
reject = array(0, dim = nTrials)
for (trial in 1:length(reject)) {  
  zk1 <- rnorm(K1);zk2 <- rnorm(K2, mean = mu, sd = sigma)  
  AUC <- Wilcoxon(zk1, zk2)  
  z <- (AUC - muAUC)/sigmaAUC
  p <- 2*pnorm(-abs(z)) # p value for individual trial
  if (p < alpha) reject[trial] = 1 
}

CI <- c(0,0); width <- -qnorm(alpha/2)
ObsvdTypeIErrRate <- sum(reject)/length(reject)
CI[1] <- ObsvdTypeIErrRate - 
  width*sqrt(ObsvdTypeIErrRate*(1-ObsvdTypeIErrRate)/nTrials)
CI[2] <- ObsvdTypeIErrRate +
 width*sqrt(ObsvdTypeIErrRate*(1-ObsvdTypeIErrRate)/nTrials)
cat("alpha = ", alpha, "\n")
cat("ObsvdTypeIErrRate = ", 
    ObsvdTypeIErrRate, "\n95% confidence interval = ", CI, "\n")
exact <- binom.test(sum(reject),n = 2000,p = alpha)
cat("exact 95% CI = ", as.numeric(exact$conf.int), "\n")
