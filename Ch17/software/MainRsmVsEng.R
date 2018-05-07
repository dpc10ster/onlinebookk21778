# mainRsmVsEng.R
rm(list = ls())
library(RJafroc)
library(ggplot2)
#library(binom)
source("FpTp2FpfTpf.R")
source("PlotBMErrBar.R")

K1 <- 500;K2 <- 700
for (Row in 1:8) {
  switch(Row,
         "1" = {seed <- 1;set.seed(seed);Lmax <- 1;mu <- 2.0;lambda <- 10;nu <- 1;zeta1 <- -1}, # Row 1
         "2" = {seed <- 1;set.seed(seed);Lmax <- 1;mu <- 2.5;lambda <- 10;nu <- 1;zeta1 <- -1}, # Row 2
         "3" = {seed <- 1;set.seed(seed);Lmax <- 1;mu <- 3.0;lambda <- 10;nu <- 1;zeta1 <- -1}, # Row 3
         "4" = {seed <- 2;set.seed(seed);Lmax <- 1;mu <- 2.5;lambda <- 10;nu <- 1;zeta1 <- -1}, # Row 4
         "5" = {seed <- 2;set.seed(seed);Lmax <- 2;mu <- 2.0;lambda <- 10;nu <- 1;zeta1 <- -1},# Row 5
         "6" = {seed <- 2;set.seed(seed);Lmax <- 2;mu <- 2.5;lambda <- 10;nu <- 1;zeta1 <- -1},# Row 6
         "7" = {seed <- 2;set.seed(seed);Lmax <- 2;mu <- 3.0;lambda <- 10;nu <- 1;zeta1 <- -1},# Row 7
         "8" = {seed <- 2;K1 <- 5000;K2 <- 7000;set.seed(seed);Lmax <- 2;mu <- 3.0;lambda <-  1;nu <- 1;zeta1 <- -1},# Row 8
         "9" = {seed <- 2;K1 <- 5000;K2 <- 7000;set.seed(seed);Lmax <- 2;mu <- 3.0;lambda <- 0.1;nu <- 1;zeta1 <- -1}# Row 9
  )
  cat("K1 = ", K1, 
      "\nK2 = ", K2, 
      "\nzeta1 = ", zeta1, 
      "\nseed = ", seed, 
      "\nLmax = ", Lmax, 
      "\nmu = ", mu, 
      "\nlambda = ", lambda, 
      "\nnu = ", nu, "\n")
  
  Lk2 <- floor(runif(K2, 1, Lmax + 1))
  nLesPerCase <- unique(Lk2);lesionDist <- array(dim = c(length(nLesPerCase), 2))
  for (i in nLesPerCase) lesionDist[i, ] <- c(i, sum(Lk2 == i)/K2)
  
  frocDataRaw  <- SimulateFrocDataset(mu, lambda, nu, I = 1, J = 1, K1, K2, lesionNum = Lk2, zeta1 = zeta1)
  rocDataRaw <- DfFroc2Roc(frocDataRaw)
  
  rocDataBinned <- DfBinDataset(rocDataRaw, desiredNumBins = 5, opChType = "ROC")
  fp <-rocDataBinned$NL[1,1,1:K1,1];tp <- rocDataBinned$LL[1,1,1:K2,1]
  ret1 <- FpTp2FpfTpf(fp, tp) 
  fpCountsOrg <- ret1$fpCounts;tpCountsOrg <- ret1$tpCounts;fpf <- ret1$fpf;tpf <- ret1$tpf;zetas <- ret1$zetas
  nBins <- length(fpCountsOrg)
  rocDataTable <- array(dim = c(2,nBins))
  rocDataTable[1,] <- fpCountsOrg;rocDataTable[2,] <- tpCountsOrg
  
  lambdaP <- UtilIntrinsic2PhysicalRSM(mu, lambda, nu)$lambdaP 
  nuP <- UtilIntrinsic2PhysicalRSM(mu, lambda, nu)$nuP
  aucs <- UtilAucsRSM(mu, lambdaP, nuP, lesionDist)
  cat("RSM-ROC-AUC = ", aucs$aucROC, "\n")
  print(rocDataTable)
  next
  # copy the last two rows of output to Eng program; delete bracket stuff leaving numbers only with spaces; select format 3
  # Run Program
  # compare to Table 1 in online appendix
  
  switch(Row, 
         # the following values were transferred from the Eng program output after analyzing data generated 
         # by mainRsmVsEng.R using the appropriate value of Row
         "1" = {a <- 1.0066;b <- 0.8182}, # Row 1
         "2" = {a <- 1.5134;b <- 0.7617}, # Row 2
         "3" = {a <- 1.9561;b <- 0.7643}, # Row 3
         "4" = {a <- 1.2324;b <- 0.7078}, # Row 4
         "5" = {a <- 1.3246;b <- 0.8715},# Row 5
         "6" = {a <- 1.5939;b <- 0.7325},# Row 6
         "7" = {a <- 2.1235;b <- 0.74086},# Row 7
         "8" = {a <- 2.4329;b <- 0.4719},# Row 8
         "9" = {a <- NA;b <- NA}# Row 9
  )
  
  print(PlotBMErrBar(a, b, rocDataTable, Row))
}