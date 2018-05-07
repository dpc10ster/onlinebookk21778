# mainEsOrdering.R
rm(list = ls())
source("AucsRsm.R")

seed <- NULL;set.seed(seed)
N <- 1000
wrong <- 0
for (k in 1:N) {
  mu1 <- runif(1,0.5,5); lambda1 <- runif(1,1,10); nu1 <- runif(1); Lmax <- 4
  mu2 <- runif(1,0.5,5); lambda2 <- runif(1,1,10); nu2 <- runif(1) 
  
  K2 <- 700;Lk2 <- floor(runif(K2, 1, Lmax + 1))
  nLesPerCase <- unique(Lk2);lesionDist <- array(dim = c(length(nLesPerCase), 2))
  for (i in nLesPerCase) lesionDist[i, ] <- c(i, sum(Lk2 == i)/K2)
  
  aucs1 <- AucsRsm(mu = mu1, lambda = lambda1, nu = nu1, lesionDist = lesionDist)
  
  aucs2 <- AucsRsm(mu = mu2, lambda = lambda2, nu = nu2, lesionDist = lesionDist)
  
  diffAFROC <- aucs1$aucAFROC - aucs2$aucAFROC
  diffROC <- aucs1$aucROC - aucs2$aucROC
  if (abs(diffAFROC) < abs(diffROC) ) {
    wrong <- wrong + 1
    }
}

cat("wrong effect size ordering fraction = ", wrong/N,"\n")

# cat("mu = ", mu, ", lambda = ", lambda, ", nu = ", nu, ", Lmax = ", 
#     Lmax, ",RSM AUC/ROC = ", aucs$aucROC, ",RSM AUC/AFROC = ", aucs$aucAFROC, "\n")
