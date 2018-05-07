# MainCalibrate.R
rm(list = ls())

library(xlsx)
library(MASS)
library(matrixcalc)
library("numDeriv")
library(RJafroc)

source("RocfitR.R")
source("ToIntegerRatings.R")
source("RunCorrocOnPairedData.R")
source("AvgParmsCovsOrg.R")
source("ParmsCovarsFittingModelOrg.R")
source("PairingType.R")
source("ReadTheData.R")

# IMPORTANT: ALWAYS chk before analysis
# IMPORTANT: ALWAYS chk before analysis
source("MySyncFile.R") 
# IMPORTANT: ALWAYS chk before analysis
# IMPORTANT: ALWAYS chk before analysis

if (Validation == TRUE) stop("Cannot be in validation mode here")
orgData <- ReadTheData(DataFile,selI,selJ)

z1ij <- orgData$NL
z2ij <- orgData$LL
I <- dim(z1ij)[1]
J <- dim(z1ij)[2]
K <- dim(z1ij)[3]
K2 <- dim(z2ij)[3]
K1 <- K - K2
z1ij <- z1ij[ , , 1:K1,]
dim(z1ij) <- c(I, J, K1)
dim(z2ij) <- c(I, J, K2)
orgData <- list(z1ij = z1ij, z2ij = z2ij)

calDir <- paste0("./REZ/", DataFile,"/","CAL/")
if (!dir.exists(calDir)) dir.create(calDir, recursive = TRUE)
fileNameCorrocOut <- paste(paste0(calDir, "CorrocOut"), "I", Istring, "J", Jstring, "BINS", DesiredNumBins, DataFile, sep = "_")
if (overWrite || !file.exists(fileNameCorrocOut)) { # as this takes a long time, let's save the results of the run
  retChiSigmachiOrg <- ParmsCovarsFittingModelOrg (DesiredNumBins, z1ij, z2ij)
  save(retChiSigmachiOrg, file = fileNameCorrocOut)
} else load(fileNameCorrocOut)

#### at this point we have a calibrated simulator; extract the values #####
doneOrg <- retChiSigmachiOrg$done
ChiOrg <- retChiSigmachiOrg$Chi
SigmaChiOrg <- retChiSigmachiOrg$SigmaChi
AxOrg <- retChiSigmachiOrg$Ax # Az arm 1
AyOrg <- retChiSigmachiOrg$Ay # Az arm 2
RhoAxAyOrg <- retChiSigmachiOrg$RhoAxAy # Rho beteeen Az's
#### END: at this point we have a calibrated simulator; extract the values #####

retAvgOrg <- AvgParmsCovsOrg(ChiOrg,SigmaChiOrg,AxOrg,AyOrg,RhoAxAyOrg,doneOrg)
ChiOrgBar <- retAvgOrg$ChiBar
SigmaChiOrgBar <- retAvgOrg$SigmaChiBar
AxOrgBar <- retAvgOrg$Ax
AyOrgBar <- retAvgOrg$Ay
RhoAxAyOrgBar <- retAvgOrg$RhoAxAy 
#### END: average the values over each of the 4 types of comparisons #####

## at this point we have a calibrated simulator

# save the values
CalSimuVals <- list(
  orgData = orgData,
  ChiOrgBar = ChiOrgBar,
  SigmaChiOrgBar = SigmaChiOrgBar,
  AxOrgBar = AxOrgBar,
  AyOrgBar = AyOrgBar,
  RhoAxAyOrgBar = RhoAxAyOrgBar
)
MyFileCalSimuVals <- paste(paste0(calDir, "CalSimuVals"), "I", Istring, "J", Jstring, "BINS", DesiredNumBins, DataFile, sep = "_")
if (overWrite || !file.exists(MyFileCalSimuVals)) {
  save(CalSimuVals, file = MyFileCalSimuVals)
}
## print out values for Excel summary
cat("ChiOrgBar = \n")
for (p in 1:6) cat(ChiOrgBar[p,],"\n") 
cat("\n")
cat("SigmaChiOrgBar = \n")
for (p in 1:6) {
  for (i1 in 1:6) {
    cat(SigmaChiOrgBar[p,i1,],"\n") # Table 
  }
  cat("\n")
}
