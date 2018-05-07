rm(list = ls())
library(mvtnorm)
source("Derivatives.R")

X <- seq(-3, 3, by = 0.1)
Y <- seq(-3, 3, by = 0.1)
MuX <- seq(-2.5, 2.5, by = 0.1)
MuY <- seq(-2.5, 2.5, by = 0.1)
Rho1 <- seq(-0.9, 0.9, by = 0.1)
Rho2 <- seq(-0.9, 0.9, by = 0.1)
Alpha <- seq(0.1, 0.9, by = 0.1)

delta <- 1e-6

dervX <- function(x, y, rho1){
  sigma <- rbind(c(1, rho1), c(rho1, 1))
  return(c(F1dx(x, y, rho1), 
           (pmvnorm(c(-Inf, -Inf), c(x + delta, y), sigma = sigma) - pmvnorm(c(-Inf, -Inf), c(x, y), sigma = sigma))/delta,
           x, y, rho1)
         )
}

derivatives1 <- data.frame(dervAnly = double(), dervNum = double(), x = double(), y = double(), rho1 = double())
for (x in X){
  for (y in Y){
    for (rho1 in Rho1){
      sigma <- rbind(c(1, rho1), c(rho1, 1))
      sigmaDelta <- rbind(c(1, rho1 + delta), c(rho1 + delta, 1))
      dervX <- c(F1dx(x, y, rho1), 
                 (pmvnorm(c(-Inf, -Inf), c(x + delta, y), sigma = sigma) - pmvnorm(c(-Inf, -Inf), c(x, y), sigma = sigma))/delta,
                 x, y, rho1)
      dervY <- c(F1dy(x, y, rho1), 
                 (pmvnorm(c(-Inf, -Inf), c(x, y + delta), sigma = sigma) - pmvnorm(c(-Inf, -Inf), c(x, y), sigma = sigma))/delta,
                 x, y, rho1)
      dervRho1 <- c(F1drho1(x, y, rho1), 
                    (pmvnorm(c(-Inf, -Inf), c(x, y), sigma = sigmaDelta) - pmvnorm(c(-Inf, -Inf), c(x, y), sigma = sigma))/delta,
                    x, y, rho1)
      tempDataFrame <- rbind(dervX, dervY, dervRho1)
      colnames(tempDataFrame) <- c("dervAnly", "dervNum", "x", "y", "rho1")
      tempDataFrame <- as.data.frame(tempDataFrame)
      derivatives1 <- rbind.data.frame(derivatives1, tempDataFrame)
    }
  }
}


derivatives2 <- data.frame(dervAnly = double(), dervNum = double(), x = double(), y = double(), rho1 = double())
for (x in X){
  for (y in Y){
    for (muX in MuX){
      for (muY in MuY){
        for (rho2 in Rho2){
          for (alpha in Alpha){
            sigma <- rbind(c(1, rho2), c(rho2, 1))
            sigmaDelta <- rbind(c(1, rho2 + delta), c(rho2 + delta, 1))
          }
        }
      }
    }
  }
}


parameters <- c(list(muIniX, muIniY, alphaIni, rhoIniNor, rhoIniAbn), as.list(zetaIniX), as.list(zetaIniY))
namesVector <- c("muX", "muY", "alpha", "rhoNor", "rhoAbn")
for (z in 1:length(zetaIniX)){
  namesVector <- c(namesVector, paste0("zetaX", z))
}
for (z in 1:length(zetaIniY)){
  namesVector <- c(namesVector, paste0("zetaY", z))
}
names(parameters) <- namesVector

CorCBMNLLNew <- AddZetaX(CorCBMNLL, length(zetaIniX))
CorCBMNLLNew <- AddZetaY(CorCBMNLLNew, length(zetaIniY))
arguments <- c(parameters, list(FPCounts = FPCounts), list(TPCounts = TPCounts))

delta <- 1e-6
arguments2 <- arguments
arguments2$zetaX1 <- arguments2$zetaX1 + delta
arguments3 <- arguments
arguments3$muX <- arguments$muX + delta
arguments4 <- arguments3
arguments4$zetaX1 <- arguments3$zetaX1 + delta


dLL1 <- (-do.call(CorCBMNLLNew, arguments2) + do.call(CorCBMNLLNew, arguments)) / (delta)
dLL2 <- (-do.call(CorCBMNLLNew, arguments4) + do.call(CorCBMNLLNew, arguments3)) / (delta)
dLLfirst <- dLL1
dLLsecond <- (dLL2 - dLL1) / delta


nBins <- dim(FPCounts)[2]
FPPrEst <- rep(NA, nBins^2)
FPComb <- FPPrEst
TPPrEst <- FPPrEst
TPComb <- FPComb
bInx <- 1
for (bX in 1:nBins){
  for (bY in 1:nBins){
    FPComb[bInx] <- FPCounts[1, bX] + FPCounts[2, bY]
    TPComb[bInx] <- TPCounts[1, bX] + TPCounts[2, bY]
    bInx <- bInx + 1
  }
}
Kdd1 <- sum(FPComb)
Kdd2 <- sum(TPComb)

zetaIniX <- c(-Inf, zetaIniX, Inf)
zetaIniY <- c(-Inf, zetaIniY, Inf)

LLdzetar <- 0
LLdmudmu <- 0
bInx <- 1

bX <- 2
for (bY in 1:nBins){
  LLdzetar <- LLdzetar + (TPCounts[1, bX] + TPCounts[2, bY]) * ((F2dx(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) - F2dx(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn)) /
                                                                  (F2(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) - F2(zetaIniX[bX], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn)
                                                                   - F2(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn) + F2(zetaIniX[bX], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn)))
  + (TPCounts[1, bX + 1] + TPCounts[2, bY]) * ((-F2dx(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) + F2dx(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn)) /
                                                 (F2(zetaIniX[bX + 2], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) - F2(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn)
                                                  - F2(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn) + F2(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn)))
  + (FPCounts[1, bX] + FPCounts[2, bY]) * ((F1dx(zetaIniX[bX + 1], zetaIniY[bY + 1], rhoIniNor) - F1dx(zetaIniX[bX + 1], zetaIniY[bY], rhoIniNor)) /
                                             (F1(zetaIniX[bX + 1], zetaIniY[bY + 1], rhoIniNor) - F1(zetaIniX[bX], zetaIniY[bY + 1], rhoIniNor)
                                              - F1(zetaIniX[bX + 1], zetaIniY[bY], rhoIniNor) + F1(zetaIniX[bX], zetaIniY[bY], rhoIniNor)))
  + (FPCounts[1, bX + 1] + FPCounts[2, bY]) * ((-F1dx(zetaIniX[bX + 1], zetaIniY[bY + 1], rhoIniNor) + F1dx(zetaIniX[bX + 1], zetaIniY[bY], rhoIniNor)) /
                                                 (F1(zetaIniX[bX + 2], zetaIniY[bY + 1], rhoIniNor) - F1(zetaIniX[bX + 1], zetaIniY[bY + 1], rhoIniNor)
                                                  - F1(zetaIniX[bX + 1], zetaIniY[bY], rhoIniNor) + F1(zetaIniX[bX + 1], zetaIniY[bY], rhoIniNor)))
}


for (bX in 1:nBins){
  for (bY in 1:nBins){
    LLdmudmu <- LLdmudmu + TPComb[bInx] * ((F2dmuXdzetar(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) - F2dmuXdzetar(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn)) /
                                             (F2(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) - F2(zetaIniX[bX], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn)
                                              - F2(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn) + F2(zetaIniX[bX], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn))
                                           - (F2dmuX(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) - F2dmuX(zetaIniX[bX], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn)
                                              - F2dmuX(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn) + F2dmuX(zetaIniX[bX], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn)) *
                                             (F2dx(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) - F2dmuXdzetar(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn)) /
                                             (F2(zetaIniX[bX + 1], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn) - F2(zetaIniX[bX], zetaIniY[bY + 1], muIniX, muIniY, alphaIni, rhoIniAbn)
                                              - F2(zetaIniX[bX + 1], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn) + F2(zetaIniX[bX], zetaIniY[bY], muIniX, muIniY, alphaIni, rhoIniAbn))^2)
    bInx <- bInx + 1
  }
}  