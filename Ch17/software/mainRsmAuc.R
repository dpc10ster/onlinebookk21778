rm(list = ls()) # mainRsmAuc.R
library(RJafroc)

seed <- 1;set.seed(seed)
mu <- 1; lambda <- 1; nu <- 1
ret <- UtilConvertIntrinsic2PhysicalRSM(
  mu, lambda, nu)
lambdaP <- ret$lambdaP;nuP <- ret$nuP
Lmax <- 1;K2 <- 700
Lk2 <- floor(runif(K2, 1, Lmax + 1))
nLesPerCase <- unique(Lk2)
lesionDist <- array(dim = c(length(nLesPerCase), 2))
for (i in nLesPerCase) 
  lesionDist[i, ] <- c(i, sum(Lk2 == i)/K2)
aucs <- UtilAucsRSM(
  mu, lambdaP = lambdaP, 
  nuP = nuP, lesionDistribution = lesionDist)
cat("mu = ", mu, 
    "\nlambda = ", lambda, 
    "\nnu = ", nu, 
    "\nLmax = ", Lmax, 
    "\nRSM AUC/ROC = ", aucs$aucROC, 
    "\nRSM AUC/AFROC = ", aucs$aucAFROC, 
    "\n")

