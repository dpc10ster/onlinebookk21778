# mainProprocDegenerate.R
rm(list = ls()) 
library(RJafroc)
library(ggplot2)
library(binom)
source("ProprocFits.R")
source("PlotRsmProp.R")
source("rsmFunctions.R")

pathName <- "PROPROCDegeneracy/"

# following are x,y coordinates of single point of degenerate data sets
fpfArr <- c(0,    0,    0.01, 0.25, 0.99)
tpfArr <- c(0.01, 0.75, 0,    1,    0)
for (i in 1:length(fpfArr)) {
  fpf <- fpfArr[i];tpf <- tpfArr[i]
  fileName <- paste0(pathName, as.character(fpf),", ",as.character(tpf),".lrc")

  rocData <- DfReadDataFile(fileName, format = "MRMC")
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
    nLesDistr <- t(rbind(as.numeric(unlist(attr(nLesDistr, "dimnames"))), 
                         as.vector(nLesDistr)))
  }
  
  retRsm <- FitRsmRoc(rocData, lesDistr = nLesDistr)
  mu <- retRsm$mu
  lambdaP <- retRsm$lambdaP
  nuP <- retRsm$nuP
  aucRsm <- retRsm$AUC
  
  fileName1 <- paste0(as.character(fpf),", ",as.character(tpf))
  retProp <- ProprocFits(fileName1)
  
  c1 <- as.numeric(retProp$c1);da <- as.numeric(retProp$da)
  area <- as.numeric(retProp$area)
  
  rho2 <- -(1-c1^2)/(1+c1^2)
  corr <- diag(2)
  corr[lower.tri(corr)] <- rho2
  corr[upper.tri(corr)] <- rho2
  lower <- rep(-Inf,2)
  upper <- c(-da/sqrt(2),0)
  mean <- rep(0,2)
  aucProproc <- pnorm(da/sqrt(2))+2*pmvnorm(lower, upper, mean, corr) 
  
  empOp <- UtilBinCountsOpPts(rocData)
  fpf <- empOp$fpf; tpf <- empOp$tpf
  compPlot <- PlotRsmProp(mu, lambdaP, nuP, nLesDistr, c1, da, fpf, tpf, 1, 1, K1, K2, 1)
  print(compPlot)
  cat("single operating point = (", fpf, tpf, ")", "\n")
  cat("PROPROC c parameter = ", c1, "\n")
  cat("PROPROC da parameter = ", da, "\n")
  cat("PROPROC fitted AUC = ", area, "\n\n")
}