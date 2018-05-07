# mainRSM.R 
# RSM fits vs. PROPROC and CBM fits; 
# Windows proproc results must be saved to MRMCRuns 
# directory prior to running this
rm(list = ls())
library(RJafroc)
library(ggplot2)
library(bbmle)
library(stats)
library(binom)

options(digits = 3)
reAnalyze <- FALSE;showPlot <- TRUE;saveProprocLrcFile <- FALSE

# included datasets
fileNames <-  c("TONY", "VD", "FR", 
                "FED", "JT", "MAG", 
                "OPT", "PEN", "NICO",
                "RUS", "DOB1", "DOB2", 
                "DOB3", "FZR")

f <- 1 # selected dataset
fileName <- fileNames[f]
# the datasets already exist as R objects
theData <- get(sprintf("dataset%02d", f))
# RSM ROC fitting needs to know lesionDistribution
lesionDistribution <- UtilLesionDistribution(theData)

rocData <- DfFroc2Roc(theData)

if (saveProprocLrcFile) {
  DfSaveDataFile(rocData, 
                 fileName = 
                   paste0(fileName,".lrc"), 
                 format = "MRMC")
}
I <- length(rocData$modalityID);J <- length(rocData$readerID)
K <- dim(rocData$NL)[3];K2 <- dim(rocData$LL)[3];K1 <- K - K2

## retrieve PROPROC parameters
csvFileName <- paste0(fileName, " proproc area pooled.csv")
sysCsvFileName <- system.file(
  paste0(
    "MRMCRuns/",fileName), csvFileName, package = "RJafroc")
if (!file.exists(sysCsvFileName)) 
  stop("Run Windows PROPROC for this dataset using VMwareFusion")
proprocRet <- read.csv(sysCsvFileName)
c1 <- matrix(data = 
               proprocRet$c, 
             nrow = length(unique(proprocRet$T)), 
             ncol = length(unique(proprocRet$R)), byrow = TRUE)
da <- matrix(data = 
               proprocRet$d_a, 
             nrow = length(unique(proprocRet$T)), 
             ncol = length(unique(proprocRet$R)), byrow = TRUE)

retFileName <- paste0("allResults", fileName) 
sysAnalFileName <- system.file(
  "ANALYZED/RSM6", retFileName, package = "RJafroc")

if (fileName %in% c("JT", "NICO", "DOB1", "DOB3")){
  binnedRocData <- DfBinDataset(rocData, desiredNumBins = 5, opChType = "ROC")
}else{
  binnedRocData <- rocData
}

if (reAnalyze || !file.exists(sysAnalFileName)){
  allResults <- list()
  AllResIndx <- 0
  for (i in 1:I){
    for (j in 1:J){
      AllResIndx <- AllResIndx + 1
      cat("f, i, j:", f, i, j, "\n")
      # fit to CBM
      retCbm <- FitCbmRoc(
        binnedRocData, trt = i, rdr = j)
      # fit to RSM, need lesionDistribution matrix
      retRsm <- FitRsmRoc(
        binnedRocData, trt = i, rdr = j, lesionDistribution)
      if (showPlot) {
        x <- allResults[[AllResIndx]]
        lesionDistribution <- x$lesionDistribution
        empOp <- UtilBinCountsOpPts(binnedRocData, trt = i, rdr = j)
        fpf <- empOp$fpf; tpf <- empOp$tpf
        compPlot <- gpfPlotRsmPropCbm(
          which(fileNames == fileName), 
          x$retRsm$mu, x$retRsm$lambdaP, x$retRsm$nuP, 
          lesionDistribution, c1[i, j], da[i, j],
          x$retCbm$mu, x$retCbm$alpha,
          fpf, tpf, i, j, K1, K2, c(1, length(fpf)))
        print(compPlot)
      }
    }
  }
  # safety comments
  # sysSavFileName <- 
  # paste0("/Users/Dev/rjafroc/inst/ANALYZED/RSM6/", retFileName)
  # save(allResults, file = sysSavFileName)
}else{
  load(sysAnalFileName)
  AllResIndx <- 0
  # cat(fileName,	"i, j, mu, lambdaP,	nuP, c,	da,	alpha,")
  # cat("muCbm,	AUC-RSM, AUC-PROPROC, AUC-CBM, chisq, p-value,  df\n")
  for (i in 1:I){
    for (j in 1:J){
      AllResIndx <- AllResIndx + 1
      x <- allResults[[AllResIndx]]
      if (showPlot) {
        empOp <- UtilBinCountsOpPts(binnedRocData, trt = i, rdr = j)
        fpf <- empOp$fpf; tpf <- empOp$tpf
        compPlot <- gpfPlotRsmPropCbm(
          which(fileNames == fileName), 
          x$retRsm$mu, x$retRsm$lambdaP, x$retRsm$nuP, 
          lesionDistribution = lesionDistribution, c1[i, j], da[i, j],
          x$retCbm$mu, x$retCbm$alpha,
          fpf, tpf, i, j, K1, K2, c(1, length(fpf)))
        compPlot <- compPlot + 
          theme(
            axis.text=element_text(size=10),
            axis.title=element_text(size=28,face="bold"))
        print(compPlot)
      }
      # follows same format as RSM6 Vs. Others.xlsx
      cat(fileName, i, j, x$retRsm$mu, x$retRsm$lambdaP, x$retRsm$nuP, 
          c1[i,j], da[i,j], 
          x$retCbm$alpha, x$retCbm$mu,
          x$retRsm$AUC, x$aucProp, x$retCbm$AUC, 
          x$retRsm$ChisqrFitStats[[1]], x$retRsm$ChisqrFitStats[[2]], 
          x$retRsm$ChisqrFitStats[[3]],"\n")
    }
  }
}