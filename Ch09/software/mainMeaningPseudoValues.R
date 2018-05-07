rm(list = ls()) # mainMeaningPseudoValues.R
library(RJafroc)
source("AvgPsdValsRatings.R")

options(digits = 3)
if (TRUE) {
  #fileName <- "Franken1.lrc"
  fileName <- "VanDyke.lrc"
  rocData <- DfReadDataFile(fileName, format = "MRMC")
} else {
  fileName <- "CXRinvisible3-20mm.xlsx"
  frocData <- DfReadDataFile(fileName, format = "JAFROC")
  rocData <- DfFrocC2HrRoc(frocData)
  rm(frocData)
}

K <- dim(rocData$NL)[3];K2 <- dim(rocData$LL)[3];K1 <- K-K2
I <- dim(rocData$NL)[1];J <- dim(rocData$NL)[2]

cat("data file = ", fileName, ", modality 1, reader 1\n")
FOM <- UtilFigureOfMerit(rocData,fom = "Wilcoxon")
cat("Number of non-diseased cases = ", K1, "\nnumber of diseased cases =  = ", K2, "\n")

for (i in 1:I) {
  for (j in 1:J) {
    rocData1T1R <- DfExtractDataset(rocData, trts=i, rdrs = j)
    zk1 <- rocData1T1R$NL[1,1,1:K1,1];zk2 <- rocData1T1R$LL[1,1,1:K2,1]
    pseudoVals <- UtilPseudoValues(rocData1T1R, fom = "Wilcoxon")
    pseudoVals <- pseudoVals[1,1,]
    psdVlsNor <- pseudoVals[1:K1]
    psdVlsAbn <- pseudoVals[(K1+1):(K1+K2)]
    if ((i == 1) & (j == 1)) {
      t1 <- table(round(psdVlsNor,3))
      cat("counts table for pseudovalues from non-diseased cases:")
      print(t1)
      cat("counts table for ratings from non-diseased cases:")
      tzk1 <- table(round(zk1,2))
      print(tzk1)
      t2 <- table(round(psdVlsAbn,3))
      cat("counts table for pseudovalues from diseased cases:")
      print(t2)
      cat("counts table for ratings from diseased cases:")
      tzk2 <- table(round(zk2,2))
      print(tzk2)
    }
    chk1 <- AvgPsdValsRatings(psdVlsNor, zk1, 1, zk2)
    chk2 <- AvgPsdValsRatings(psdVlsAbn, zk2, 2, zk1)
    if (abs(chk1$avgY - chk2$avgY) > 1e-6) stop("error 1")
    if (abs(chk2$avgY - FOM[i,j]) > 1e-6) stop("error 2")
    ##cat("sel. Trt. = ", i, ", sel. Rdr. = ", j, ", avg psedovalue = ", chk1$avgY, ", avg rating = ", chk1$avgR, "\n")
    next
  }
}