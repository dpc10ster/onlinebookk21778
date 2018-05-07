# mainRsmVsRocfitR.R
rm(list = ls())
library(RJafroc)
library(ggplot2)
library("numDeriv")
source("Transforms.R")
source("LL.R")
source("RocOperatingPointsFromRatingsTable.R")
source("VarianceAz.R");
source("ChisqrGoodnessOfFit.R")
source("RocfitR.R")
source("AucsRsm.R")
source("PlotRSMBM.R")

lambda <- 1;zeta1 <- -1;nBins <- 5;K1 <- 500;K2 <- 700
cat("K1 = ", K1, ", K2 = ", K2, "\n")
muArr <- c(2,3);nuArr <- c(0.15,0.25);
cat("lambda = ", lambda, ", zeta1 = ", zeta1, "\n")
seedArr <- c(2,3)
for (s in 1:2){
  seed <- seedArr[s];set.seed(seed)
  cat("seed = ", seed, "\n")
  Lmax <- 1;Lk2 <- floor(runif(K2, 1, Lmax + 1))
  nLesPerCase <- unique(Lk2)
  lesionDist <- array(dim = c(length(nLesPerCase), 2))
  for (i in nLesPerCase) lesionDist[i, ] <- c(i, sum(Lk2 == i)/K2)
  cat("Lmax = ", Lmax, "\n")
  
  for (i in 1:2){
    for (j in 1:2){
      mu <- muArr[i];nu <- nuArr[j]
      RowString <- toString(c(seed,mu,nu))
       
      while(1){      
        frocDataRaw  <- SimulateFrocDataset (
          mu, lambda, nu, 
          I = 2, J = 2, K1 = K1, K2, 
          lesionNum = Lk2, zeta1 = zeta1)
        rocDataRaw <- DfFroc2Roc(frocDataRaw)
        
        rocDataBinned <- DfBinDataset(
          rocDataRaw,
          desiredNumBins = nBins, 
          opChType = "ROC")
        if (length(unique(rocDataBinned$LL[i,j,,1])) != nBins) next
        if (length(unique(rocDataBinned$NL[i,j,1:K1,1])) != nBins) next
        break
      }      
      
      rocDataTable <- array(dim = c(2,nBins))
      rocDataTable[1,] <- table(rocDataBinned$NL[i,j,1:K1,1])
      rocDataTable[2,] <- table(rocDataBinned$LL[i,j,1:K2,1])
      
      retRocfit <- RocfitR(rocDataTable)
      if (length(retRocfit) != 6) stop("rocfit failed")
      
      aucs <- AucsRsm(mu = mu, 
                      lambda = lambda, 
                      nu = nu, 
                      lesionDist = lesionDist)
      cat("mu=", mu, 
          ",nu=", nu, 
          ",RSM-ROC-AUC = ", aucs$aucROC, 
          ",Az=", retRocfit$Az, 
          ",a=", retRocfit$a, 
          ",b =", retRocfit$b, "\n")
      compPlot <- PlotRSMBM(
        retRocfit$a, 
        retRocfit$b, 
        mu, 
        lambda, 
        nu,  lesionDist, RowString)
      print(compPlot$compPlot)
      next
    }
  }
}
