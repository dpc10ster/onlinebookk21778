rm(list = ls()) #mainCorCBM.R
library(ggplot2)
library(RJafroc)
library(mvtnorm)
library(bbmle)
library(grid)
library(binom)
library(numDeriv)

source("CorCBM.R")
source("Transforms.R")
source("addArguments.R")
source("RunCorrocOnPairedData.R")
source("PlotCORCBMCORROC.R")
source("EmpEstProbCORCBM.R")
source("EmpEstProbCORROC.R")

dataFile <- "JT"
if (dataFile == "VD"){
  orgData <- ReadDataFile("./datasets/VD.lrc", format = "MRMC")
}else if (dataFile == "FR"){
  orgData <- ReadDataFile("./datasets/FR.lrc", format = "MRMC")
}else if (dataFile == "JT"){
  orgData <- ReadDataFile("./datasets/JT.xlsx"); orgData <- FROC2HrROC(orgData)
  orgData <- ExtractDataset(orgData, 1:2, 1); orgData <- BinDataset(orgData, 5)
}
cat("Data file is", dataFile,"\n")

# seed <- 1
# set.seed(seed)
# z1z2 <- BvBnSimulator(1.250567, 1.890867, 0.3451667, 0.6144667, 0.2536000, 0.7562000, 50, 50 )
# z1 <- z1z2$z1
# z2 <- z1z2$z2
# z1 <- array(t(z1), dim = c(1, 2, 50))
# z2 <- array(t(z2), dim = c(1, 2, 50))
# simuData <- ToRJafrocDataset(z1, z2, "ROC")
# simuDataBin <- BinDataset(simuData, 5)
# orgData <- simuDataBin$dataset

aucArray <- FigureOfMerit(orgData, fom = "Wilcoxon")
I <- dim(orgData$NL)[1]
J <- dim(orgData$NL)[2]
K <- dim(orgData$NL)[3]
K2 <- dim(orgData$LL)[3]
K1 <- K - K2
FP <- orgData$NL[,,1:K1,1]
TP <- orgData$LL[,,,1]
dim(FP) <- c(I, J, K1)
dim(TP) <- c(I, J, K2)

scores <- unique(c(FP, TP)) # dpc possible problem if no data in one bin??
nBins <- length(scores)

retCORCBM <- NULL
retCORROC <- NULL
retFile <- paste0("retCORCBMCORROC_", dataFile)
if (file.exists(retFile)){
  load(retFile)
}

FPF <- array(dim = c(I, J, nBins))
TPF <- array(dim = c(I, J, nBins))
for (j in 1:J){
  iX <- 1
  jX <- j
  iY <- 2
  jY <- j
  
  muMin <<- 0
  muMax <<- 3
  alphaMin <<- 0
  alphaMax <<- 1
  rhoMin <<- -1
  rhoMax <<- 1
  zetaMin <<- -3
  zetaMax <<- 5.5
  
  aucX <- aucArray[iX, jX]
  aucY <- aucArray[iY, jY]
  
  z1X <- FP[iX, jX, ]
  z1Y <- FP[iY, jY, ]
  FPCounts <- array(dim = c(nBins, nBins))
  binnedFPX <- as.integer(cut(z1X, c(scores, Inf), right = FALSE)) # dpc why is this needed? seems to have no effect
  binnedFPY <- as.integer(cut(z1Y, c(scores, Inf), right = FALSE))
  
  z2X <- TP[iX, jX, ]
  z2Y <- TP[iY, jY, ]
  TPCounts <- array(dim = c(nBins, nBins))
  binnedTPX <- as.integer(cut(z2X, c(scores, Inf), right = FALSE))
  binnedTPY <- as.integer(cut(z2Y, c(scores, Inf), right = FALSE))
  
  for (bX in 1:nBins){
    for (bY in 1:nBins){
      FPCounts[bX, bY] <- sum((binnedFPX == bX & binnedFPY == bY))
      TPCounts[bX, bY] <- sum((binnedTPX == bX & binnedTPY == bY))
    }
  }
  rhoIniNor <- cor(z1X, z1Y)
  rhoIniAbn <- cor(z2X, z2Y)
  
  fpX <- rowSums(FPCounts)
  tpX <- rowSums(TPCounts)
  fpY <- colSums(FPCounts)
  tpY <- colSums(TPCounts)
  FPF[1, j, ] <- cumsum(rev(fpX))/sum(fpX)
  TPF[1, j, ] <- cumsum(rev(tpX))/sum(tpX)
  
  FPF[2, j, ] <- cumsum(rev(fpY))/sum(fpY)
  TPF[2, j, ] <- cumsum(rev(tpY))/sum(tpY)
  
  if (!file.exists(retFile)){
    retCORCBM <- c(retCORCBM, list(CorCBM(FPCounts, TPCounts, aucX, aucY, K1, K2, rhoIniNor, rhoIniAbn)))
    
    z1b <- rbind(z1X, z1Y)
    z2b <- rbind(z2X, z2Y)
    retCORROC <- c(retCORROC, list(RunCorrocOnPairedData(z1b, z2b, 5)))
  }
}

empPrCORCBM <- NULL
estPrCORCBM <- NULL
empPrCORROC <- NULL
estPrCORROC <- NULL
pValCORCBM <- NULL
pValCORROC <- NULL

if (dataFile == "VD"){
  ciPoint <- rbind(c(3, 2),
                   c(3, 3),
                   c(2, 3),
                   c(3, 3),
                   c(3, 3))
}else if (dataFile =="FR"){
  ciPoint <- rbind(c(2, 2),
                   c(3, 3),
                   c(3, 3),
                   c(2, 2))
}

for (j in 1:J){
  fpf <- FPF[ , j, ]
  tpf <- TPF[ , j, ]
  p <- PlotCORCBMCORROC(retCORCBM[[j]], retCORROC[[j]], fpf, tpf, ciPoint[j, ])
  print(p[[1]])
  p1 <- p[[1]]
  gt <- ggplot_gtable(ggplot_build(p1))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid.draw(gt)
  
  print(p[[2]])
  p2 <- p[[2]]
  gt <- ggplot_gtable(ggplot_build(p2))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid.draw(gt)
}

options(digits = 4) # output significant digits
cat("\n")
cat("Result table for CORCBM:")
cat("\n\n")
for (j in 1:J){
  cat("j =", j, "\n")
  cat("Parameters:", retCORCBM[[j]]$muX, retCORCBM[[j]]$alphaX, retCORCBM[[j]]$muY, retCORCBM[[j]]$alphaY, retCORCBM[[j]]$rhoNor, retCORCBM[[j]]$rhoAbn2, retCORCBM[[j]]$aucX, retCORCBM[[j]]$aucY, "\n")
  cat("StdErrs:", retCORCBM[[j]]$stdErr, "\n")
  cat("Area test statistic:", retCORCBM[[j]]$aStat, "\n")
  cat("P-Value:", retCORCBM[[j]]$aPval)
  cat("\n\n")
}

cat("\n")
cat("Result table for CORROC:")
cat("\n\n")
for (j in 1:J){
  cat("j =", j, "\n")
  cat("Parameters:", c(retCORROC[[j]]$parms, retCORROC[[j]]$Ax, retCORROC[[j]]$Ay), "\n")
  vars <- diag(retCORROC[[j]]$Cov)
  stdErrCORROC <- sqrt(vars)
  aX <- retCORROC[[j]]$parms[1]
  bX <- retCORROC[[j]]$parms[2]
  daX <- dnorm(aX / sqrt(1 + bX^2)) / sqrt(1 + bX^2)
  dbX <- -aX * bX * dnorm(aX / sqrt(1 + bX^2)) / (1 + bX^2)^(3/2)
  stdErrAx <- sqrt(daX^2 * vars[1] + dbX^2 * vars[2] + daX * dbX * retCORROC[[j]]$Cov[1, 2])
  aY <- retCORROC[[j]]$parms[3]
  bY <- retCORROC[[j]]$parms[4]
  daY <- dnorm(aY / sqrt(1 + bY^2)) / sqrt(1 + bY^2)
  dbY <- -aY * bY * dnorm(aY / sqrt(1 + bY^2)) / (1 + bY^2)^(3/2)
  stdErrAy <- sqrt(daY^2 * vars[3] + dbY^2 * vars[4] + daY * dbY * retCORROC[[j]]$Cov[3, 4])
  stdErrCORROC <- c(stdErrCORROC, stdErrAx, stdErrAy)
  cat("StdErrs:", stdErrCORROC, "\n")
  
  derivs <- c(daX, dbX, -daY, -dbY)
  varDiff <- 0
  for (k in 1:4){
    for (l in 1:4){
      varDiff <- derivs[l] * derivs[k] * retCORROC[[j]]$Cov[k, l] + varDiff
    }
  }
  aStat <- abs(retCORROC[[j]]$Ax - retCORROC[[j]]$Ay)/sqrt(varDiff)
  aPval <- 1 - pnorm(aStat)
  cat("Area test statistic:", aStat, "\n")
  cat("P-Value:", aPval)
  cat("\n\n")
}

