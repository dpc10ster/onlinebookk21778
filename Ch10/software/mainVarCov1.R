#mainVarCov1.R
rm(list = ls())
library(RJafroc)
source("VarCov1Bs.R")
source("Wilcoxon.R")
source("VarCov1Bs.R")
source("VarCov1Jk.R")
source("VarCovs.R")
source("VarCovMtrxDLStr.R")

seed <- 1;set.seed(seed)
fileName <- "CXRinvisible3-20mm.xlsx"

frocData <- DfReadDataFile(fileName, format = "JAFROC", newExcelFileFormat = FALSE)
rocData <- DfFroc2Roc(frocData)
jSelect <- 1  # selects the reader to be analyzed
rocData1R <- DfExtractDataset(rocData, rdrs = jSelect)

zik1 <- rocData1R$NL[,1,,1];K <- dim(zik1)[2]
I <- dim(zik1)[1]
zik2 <- rocData1R$LL[,1,,1];K2 <- dim(zik2)[2]
K1 <- K-K2;zik1 <- zik1[,1:K1]

FOM <- array(dim=c(I))
for (i in 1:I) {
  FOM[i] <- Wilcoxon(zik1[i,],zik2[i,])    
}

FOMik <- array(dim = c(I,  K))
for (i in 1:I) {
  for (k in 1:K) {
    if (k <= K1) {
      FOMik[i,k] <- Wilcoxon(zik1[i,-k], zik2[i,])
    } else {
      FOMik[i,k] <- 
        Wilcoxon(zik1[i,], zik2[i,-(k-K1)])
    } 
  }
}
ret1 <- VarCov1Jk(FOMik)

cat("data file = ", fileName, "\n")
cat("number of treatments = ", I, 
    "\nnumber of non-diseased cases = ", K1, 
    "\nnumber of diseased cases = ", K2, "\n")
cat("reader = ", jSelect, "\n")
cat("OR variance components using jackknife\n")
cat("Variance = ",  ret1$Var, 
    "\nCov1 = ",  ret1$Cov1, 
    "\nrho = ",  ret1$Cov1/ret1$Var, "\n")
#to save the bs Auc values
B <- 2000;aucBs <- array(dim = c(I,B))
for (b in 1 : B){
  # bs indices for non-diseased    
  k1b <- ceiling( runif(K1) * K1 )
  # bs indices for diseased  
  k2b <- ceiling( runif(K2) * K2 )
  for ( i in 1 : I) aucBs[i,b] <- 
    Wilcoxon(zik1[i,k1b], zik2[i,k2b]) 
}

ret2 <- VarCov1Bs(aucBs)
cat("OR variance components using bootstrap\n")
cat("Variance = ",  ret2$Var, "\nCov1 = ",  
    ret2$Cov1, "\nrho = ",  
    ret2$Cov1/ret2$Var, "\n")

mtrxDLStr <- VarCovMtrxDLStr(rocData1R)
VarCovDLStr <- VarCovs(mtrxDLStr)

cat("OR variance components using DeLong method\n")
cat("Variance = ",  VarCovDLStr$var, 
    "\nCov1 = ",  
    VarCovDLStr$cov1, "\nrho = ",  
    VarCovDLStr$cov1/VarCovDLStr$var, "\n")

