rm(list = ls()) # mainStatPower.R
source("Wilcoxon.R");source("EffectSize.R")

seed <- 1;set.seed(seed)
mu <- 1.5;muAH <- 2.1;sigma <- 1.3
K1 <- 50;K2 <- 52#;K1 <- K1*2;K2 <- K2*2

# cheat to find the population mean and std. dev.
AUC <- array(dim = 10000)
for (i in 1:length(AUC)) {
  zk1 <- rnorm(K1)
  zk2 <- rnorm(K2, mean = mu, sd = sigma)  
  AUC[i] <- Wilcoxon(zk1, zk2)
}
sigmaAUC  <-  sqrt(var(AUC));muAUC   <-  mean(AUC)

T <- 2000
alpha <- 0.05 # size of test
reject = array(0, dim = T)
for (t in 1:length(reject)) {  
  zk1 <- rnorm(K1)
  zk2 <- rnorm(K2, mean = muAH, sd = sigma)  
  AUC <- Wilcoxon(zk1, zk2)  
  obsvdZ <- (AUC - muAUC)/sigmaAUC
  p <- 2*pnorm(-abs(obsvdZ)) # p value for individual t
  if (p < alpha) reject[t] = 1 
}
 
ObsvdTypeIErrRate <- sum(reject)/length(reject)
CI <- c(0,0);width <- -qnorm(alpha/2)
CI[1] <- ObsvdTypeIErrRate - 
  width*sqrt(ObsvdTypeIErrRate*(1-ObsvdTypeIErrRate)/T)
CI[2] <- ObsvdTypeIErrRate + 
  width*sqrt(ObsvdTypeIErrRate*(1-ObsvdTypeIErrRate)/T)
cat("alpha = ", alpha, "\n")
cat("#non-diseased images = ", K1, 
    "\n#diseased images = ", K2, "\n")
cat("obsvdPower = ", ObsvdTypeIErrRate, "\n")
cat("95% confidence interval = ", CI, "\n")
cat("Effect Size = ", EffectSize(mu, sigma, muAH, sigma), "\n")
