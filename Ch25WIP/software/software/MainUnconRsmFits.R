rm(list = ls()) # MainUnconRsmFits.R redid tony analysis, dec 14 2016
library(RJafroc)
library(bbmle)
library(stats)
library(ggplot2)
require("mvtnorm")

source("OperatingCharacteristics.R")
source("Transforms.R")
source("optFunctions.R")
source("addArguments.R")
source("FitOldROCCurve.R")
source("PlotRsmProp.R")

UNINITIALIZED <- -Inf
minZeta <<- -20
maxZeta <<- 20
minLambdaP <<- 0.001
maxLambdaP <<- 10
minNuP <<- 0
maxNuP <<- 1
minMu <<- 0.001
maxMu <<- 20 

frocData <- ReadDataFile("ALMLCWBZSZ_Finale_20100402.xlsx")  
rocData <- FROC2HrROC(frocData)

lesionNum <- frocData$lesionNum
nLesDistr <- table(lesionNum)
if (length(nLesDistr) == 1) {
  nLesDistr <- c(lesionNum[1], 1)
  dim(nLesDistr) <- c(1, 2)
}else{
  nLesDistr <- t(rbind(as.numeric(unlist(attr(nLesDistr, "dimnames"))), as.vector(nLesDistr)))
}
I <- length(rocData$modalityID)
J <- length(rocData$readerID)

# following values for proproc parameters are from file 
c <- rbind(c(-0.132280361, -0.086965135, -0.1444418524, 0.080460161, 0.222558765),
           c(-0.081742476, 0.049764485, -0.132612623, 0.118222633, 0.078103299))
da <- rbind(c(1.197239295, 1.771175637, 1.481934876, 1.513756914, 1.740157217),
            c(0.628125133, 0.973878556, 1.155870661, 1.620175716, 0.89288159))

proprocAuc <- c*0

for (i in 1:I){
  for (j in 1:J){
    rho2 <- -(1-c[i,j]^2)/(1+c[i,j]^2)
    corr <- diag(2)
    corr[lower.tri(corr)] <- rho2
    corr[upper.tri(corr)] <- rho2
    lower <- rep(-Inf,2)
    upper <- c(-da[i,j]/sqrt(2),0)
    mean <- rep(0,2)
    proprocAuc[i,j] <- pnorm(da[i,j]/sqrt(2))+2*pmvnorm(lower, upper, mean, corr)
    rsmRet <- FitOldROCCurve(rocData, i, j, nLesDistr)
    empOp <- EmpiricalOpCharac(rocData, i, j, opChType = "ROC")$ROCPoints
    mu <- as.numeric(rsmRet$mu$mu)
    lambdaP <- as.numeric(rsmRet$lambdaP$lambdaP)
    nuP <- as.numeric(rsmRet$nuP$nuP)
    aucSM <- PlotRsmOperatingCharacteristics(mu = mu, lambda = lambdaP * mu, nu = -log(1 - nuP)/mu, type = "ROC")$aucROC
    fpf <- empOp$FPF; tpf <- empOp$TPF
    fpf <- fpf[-c(1, length(fpf))]
    tpf <- tpf[-c(1, length(tpf))]
    compPlot <- PlotRsmProp(mu, lambdaP, nuP, nLesDistr, c[i, j], da[i, j], fpf, tpf, i, j)
    
    #print(compPlot)
    cat("I =", i, ", J =", j, ", lambdaP =", lambdaP, ", nuP =", nuP, ", mu", mu,", AUCSm =", aucSM, 
        ", proprocAuc =", proprocAuc[i,j], "\n")
  }
}
