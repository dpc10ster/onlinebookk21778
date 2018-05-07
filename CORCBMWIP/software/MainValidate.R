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

# IMPORTANT: ALWAYS chk before analysis
# IMPORTANT: ALWAYS chk before analysis
source("MySyncFile.R") 
# IMPORTANT: ALWAYS chk before analysis
# IMPORTANT: ALWAYS chk before analysis

if (Validation == FALSE) stop("Cannot be in calibration mode here")

CalibrationFile <- paste(paste0("./REZ/",DataFile,"/","CAL/","CalSimuVals"), "I", Istring, "J", Jstring, "BINS", 
                         DesiredNumBins, DataFile, sep = "_")
load(CalibrationFile)
ChiOrgBar <- CalSimuVals$ChiOrgBar
SigmaChiOrgBar <- CalSimuVals$SigmaChiOrgBar
z1ij <- CalSimuVals$orgData$z1ij;z2ij <- CalSimuVals$orgData$z2ij
calI <- dim(z1ij)[1];calJ <- dim(z1ij)[2];calK1 <- dim(z1ij)[3];calK2 <- dim(z1ij)[4]
if (nhCondition == TRUE) {
  nh <- CreateNHParameters(calI, calJ, ChiOrgBar, SigmaChiOrgBar)
  ChiOrgBar <- nh$ChiOrgBar
  SigmaChiOrgBar <- nh$SigmaChiOrgBar
}
## at this point we have a calibrated simulator

ret <- BootstrapCIs (nBoots, z1ij, z2ij, DesiredNumBins, c(1:calI), aucEst)
CI <- ret$CI; bsMean <- ret$bsMean; bsVar <- ret$bsVar
## at this point we have bs CIs and Variance

## change numbers of readers and cases for simulations
if (!is.null(simuI)){
  I <- simuI
} else I <- calI
if (!is.null(simuJ)){
  J <- simuJ
} else J <- calJ
if (!is.null(simuK1)){
  K1 <- simuK1
} else K1 <- calK1
if (!is.null(simuK2)){
  K2 <- simuK2
} else K2 <- calK2
## change numbers of readers and cases for simulations

if (simSelect != "RM") {
  retSimulator <- setupSimulator (S, I, J, K1, K2, ChiOrgBar, SigmaChiOrgBar, orgData, CalSimuVals, 
                                nhCondition, resampleChi, DesiredNumBins, RMVarStr, aucEst = aucEst) # TEMP*DPC
  ChiNew <- retSimulator$ChiNew
} else {
  retSimulator <- setupSimulator (S, I, J, K1, K2, ChiOrgBar, SigmaChiOrgBar, orgData, CalSimuVals,
                                  nhCondition, resampleChi, DesiredNumBins, RMVarStr, aucEst = aucEst) # TEMP*DPC
  ChiNew <- NA
}
if (simSelect == "RM") 
{
  varComp <- retSimulator$varComp
  mu <- retSimulator$mu
  tau <- retSimulator$tau
}
MyFile <- retSimulator$MyFile
set.seed(seed)
saveSim <- simulateTheData (simSelect, MyFile, S, I, J, K1, K2, ChiOrg, ChiOrgBar, 
                            mu, tau, ChiNew, varComp, DesiredNumBins, overWrite, isBinned = TRUE)
## at this point we have the simulated data; ready to plot histograms and values

PlotsValues(I, saveSim, CI, bsMean, bsVar, nhCondition)
