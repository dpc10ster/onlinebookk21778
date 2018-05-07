# RSM fits vs. PROPROC and CBM fits; #LOCK line numbers from now on 
# Windows proproc results must be saved to MRMCRuns directory prior to running this
rm(list = ls()) # mainRSM.R 
library(RJafroc)
library(bbmle)
library(stats)
library(ggplot2)
library(binom)
library("caTools")
require("mvtnorm")

source("lesionDistribution.R")
source("rsmFunctions.R")
source("PlotRsmPropCbm.R")
source("loadDataFile.R")
source("ProprocFits.R") # contains results of PROPROC fits run on Windows machine

pathName <- "../../../06 E Online Appendices/E24 Datasets/"
# included datasets
fileNames <-  c("TONY", "VD", "FR", "FED", "JT", "MAG", "OPT", "PEN", "NICO",
                "RUS", "DOB1", "DOB2", "DOB3", "FZR")
pVal <- NULL
for (f in 1:length(fileNames)){
  cat("fileName = ", fileNames[f],"\n")
  retFileName <- paste0("ANALYZED/RSM-PROPROC-CBM/", "saveRetRoc", fileNames[f])
  if (file.exists(retFileName)){
    load(retFileName)
    load(retFileName)
    frocData <- loadDataFile(fileNames[f], pathName)
    I <- length(frocData$modalityID)
    J <- length(frocData$readerID)
    s <- 1
    for (i in 1:I){
      for (j in 1:J){
        pVal <- c(pVal, retSmRoc[[s]]$gdnss$gdnss)
        s <- s + 1
      }
    } 
  } else {
    stop("missing file")
  }
}
naCount <- length(pVal[is.na(pVal)])
pVal <- pVal[!is.na(pVal)]

df <- data.frame(value = pVal)
histogram <- ggplot(df, aes(x = value)) + geom_histogram(color = "white", binwidth = 0.1) + xlab("pValue")
print(histogram)
