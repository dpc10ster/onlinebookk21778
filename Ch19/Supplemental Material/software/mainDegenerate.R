rm(list = ls()) # mainDegenerate.R
library(mvtnorm)
library(RJafroc)
library(ggplot2)
library(binom)
source("ProprocFits.R")
source("PlotRsmProp.R")
source("rsmFunctions.R")
# following are x,y coordinats of single point of degenerate data
fileNames <- c("0, 0.01", "0, 0.75", "0.01, 0", "0.25, 1", "0.99, 0")
fileName <- fileNames[5]
filePath <- paste0("./MRMCRuns/", fileName, ".lrc")

rocData <- DfReadDataFile(filePath, format = "MRMC")
I <- length(rocData$modalityID)
J <- length(rocData$readerID)
K <- dim(rocData$NL)[3]
K2 <- dim(rocData$LL)[3]
K1 <- K - K2
lesionNum <- rocData$lesionNum
nLesDistr <- table(lesionNum)
if (length(nLesDistr) == 1) {
  nLesDistr <- c(lesionNum[1], 1)
  dim(nLesDistr) <- c(1, 2)
}else{
  nLesDistr <- t(rbind(as.numeric(unlist(attr(nLesDistr, "dimnames"))), as.vector(nLesDistr)))
}


retRsm <- FitRsmRoc(rocData, 1, 1)
mu <- retRsm$mu$mu
lambdaP <- retRsm$lambdaP$lambdaP
nuP <- retRsm$nuP$nuP
aucRsm <- retRsm$AUC$AUC

retProp <- ProprocFits(fileName)

c1 <- retProp$c1;da <- retProp$da

rho2 <- -(1-c1^2)/(1+c1^2)
corr <- diag(2)
corr[lower.tri(corr)] <- rho2
corr[upper.tri(corr)] <- rho2
lower <- rep(-Inf,2)
upper <- c(-da/sqrt(2),0)
mean <- rep(0,2)
aucProproc <- pnorm(da/sqrt(2))+2*pmvnorm(lower, upper, mean, corr) 

empOp <- GetOperatingPoints(rocData, 1, 1, opChType = "ROC")
fpf <- empOp$FPF; tpf <- empOp$TPF
compPlot <- PlotRsmProp(mu, lambdaP, nuP, nLesDistr, c1, da, fpf, tpf, 1, 1, K1, K2, 1)
print(compPlot)
