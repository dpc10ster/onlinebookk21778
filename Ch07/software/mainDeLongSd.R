rm(list = ls()) # mainDeLongSd.R
source("Wilcoxon.R");source("DeLongVar.R")

seed <- 1;set.seed(seed)
mu <- 1.5;sigma <- 1.3;K1 <- 50;K2 <- 52
cat("seed = ", seed, "\nK1 = ", K1, "\nK2 = ", K2, 
    "\nmu = ", mu, "\nsigma = ", sigma, "\n")

# brute force method to find the population mean and stdDev
empAuc <- array(dim = 10000)
for (i in 1:length(empAuc)) {
  zk1 <- rnorm(K1);zk2 <- rnorm(K2, mean = mu, sd = sigma)  
  empAuc[i] <- Wilcoxon(zk1, zk2)
}
meanempAuc   <-  mean(empAuc)
stdDevempAuc  <-  sqrt(var(empAuc))
cat("population mean empAuc = ", meanempAuc, 
    "\npopulation stdDev empAuc = ", stdDevempAuc, "\n")

# one more trial
zk1 <- rnorm(K1);zk2 <- rnorm(K2, mean = mu, sd = sigma)  
empAuc <- Wilcoxon(zk1, zk2)
ret  <- DeLongVar(zk1,zk2)
stdDevDeLong <- sqrt(ret)
cat("1 sample empAuc = ", empAuc, 
    "\nstdDev DeLong = ", stdDevDeLong, "\n")
