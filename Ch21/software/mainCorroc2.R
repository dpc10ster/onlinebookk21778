rm(list = ls()) # mainCorroc2.R
library(RJafroc)
source("RunCorrocOnPairedData.R")

jSelected <- 1
dataFile <- "VanDyke.lrc"
dataset <- DfReadDataFile(dataFile, format = "MRMC")
dataset <- DfExtractDataset(dataset, rdrs=jSelected)
K1 <- dim(dataset$NL)[3]-dim(dataset$LL)[3]
K2 <- dim(dataset$LL)[3]
z1b <- dataset$NL[,1,1:K1,1]
z2b <- dataset$LL[,1,1:K2,1]
retCorrocii <- RunCorrocOnPairedData (z1b,z2b,5)
if (length(retCorrocii) == 1) stop("CorrocII failed") else # skip degenerate readers
{ 
  parametersOrig <- retCorrocii$parms
  Cov <- retCorrocii$Cov
  cat("The 6 parameters are ", parametersOrig, "\n")
  cat("The 2 sided pValue is ", retCorrocii$pValue, "\n")
  cat("The covariance matrix is:\n")
  print(Cov)
}
