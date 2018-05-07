rm(list = ls()) #mainHT1R1M.R
source("Wilcoxon.R")

seed <- 1;set.seed(seed)
mu <- 1.5;sigma <- 1.3;K1 <- 50;K2 <- 52

# cheat to find the population mean and std. dev.
AUC <- array(dim = 10000)
for (i in 1:length(AUC)) {
  zk1 <- rnorm(K1)
  zk2 <- rnorm(K2, mean = mu, sd = sigma)  
  AUC[i] <- Wilcoxon(zk1, zk2)
}
meanAUC   <-  mean(AUC);sigmaAUC  <-  sd(AUC)
cat("pop mean AUC = ", meanAUC, 
    "\npop sigma AUC = ", sigmaAUC, "\n")

# one more trial, this is the one we want 
# to compare to meanAUC, 
zk1 <- rnorm(K1);zk2 <- rnorm(K2, mean = mu, sd = sigma) 
AUC <- Wilcoxon(zk1, zk2)
cat("New AUC = ", AUC, "\n")

z <- (AUC - meanAUC)/sigmaAUC
#z <- qnorm(0.05/2)
cat("z-statistic = ", z, "\n")

# p value for two-sided AH
p2tailed <- pnorm(-abs(z)) + (1-pnorm(abs(z)))
# p value for one-sided AH > 0
p1tailedGT <- 1-pnorm(z)
# p value for one-sided AH < 0 
p1tailedLT <- pnorm(z)
alpha  <- 0.05

# critical value for two-sided AH: 
# AUC not equal to meanAUC
z2tailed <- -qnorm(alpha/2)
# critical value for one-sided AH: 
# AUC > meanAUC
z1tailedGT <- qnorm(1-alpha)
# critical value for one-sided AH: 
# AUC < meanAUC
z1tailedLT <- qnorm(alpha)

cat("alpha of test = ", alpha, "\n")
cat("\nTwo-sided AH: AUC not equal to meanAUC", "\n")
cat("Critical value for two-sided AH:", z2tailed, "\n")
cat("p value for two-sided AH:", p2tailed, "\n")

cat("\nOne-sided AH: AUC > meanAUC", "\n")
cat("Critical value for one-sided AH:", z1tailedGT, "\n")
cat("p value for two-sided AH:", p1tailedGT, "\n")

cat("\nOne-sided AH: AUC < meanAUC", "\n")
cat("Critical value for one-sided AH:", z1tailedLT, "\n")
cat("p value for two-sided AH:", p1tailedLT, "\n")
