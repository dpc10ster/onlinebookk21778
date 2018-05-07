# MainValidate.R 
rm(list = ls())

library(xlsx)
library(MASS)
library(lmf)
library(matrixcalc)
library("numDeriv")
library(foreach)
library(doParallel)
library(doMC)
library(doRNG)
library(ggplot2)
library(RJafroc)

source("ToIntegerRatings.R")
source("orCovariances.R")
source("RunCorrocOnPairedData.R")
source("SimulateRocDataNew.R")
source("BootstrapCIs.R")
source("PairingType.R")
source("RocfitR.R")
source("RoeMetzVarStr.R")
source("RMSimulator.R")
source("SigmaChiNew.R")
source("PlotsValues.R")
source("setupSimulator.R")
source("simulateTheData.R")
source("CreateNHParameters.R")
source("SelectRMVarStr.R")
source("Cov2Validation.R")

# IMPORTANT: ALWAYS chk before analysis
# IMPORTANT: ALWAYS chk before analysis
source("MySyncFile.R") 
# IMPORTANT: ALWAYS chk before analysis
# IMPORTANT: ALWAYS chk before analysis

CalibrationFile <- paste(paste0("./REZ/",DataFile,"/","CAL/","CalSimuVals"), "I", Istring, "J", Jstring, "BINS", 
                         DesiredNumBins, DataFile, sep = "_")
load(CalibrationFile)
z1ij <- CalSimuVals$orgData$z1ij;z2ij <- CalSimuVals$orgData$z2ij
calI <- dim(z1ij)[1];calJ <- dim(z1ij)[2];calK1 <- dim(z1ij)[3];calK2 <- dim(z1ij)[4]

cov2 <- Cov2Validation(nBoots, z1ij, z2ij, DesiredNumBins, c(1:calI), aucEst)

if (!is.null(simuJ)){
  J <- simuJ
} else J <- calJ
cov2[cov2[, (J + 1)] > 1.5 *  cov2[, (J + 2)],1:J]

cov2[cov2[, (J + 1)] <= 1.5 *  cov2[, (J + 2)],1:J]
