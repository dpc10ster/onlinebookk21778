rm(list = ls()) #mainOrhDbmh.R
library(RJafroc)
ROC <- FALSE
if (ROC) {
  #fileName <- "Franken1.lrc"
  fileName <- "VanDyke.lrc"
  rocData <- DfReadDataFile(
    fileName, 
    format = "MRMC")
} else {
  fileName <- "CXRinvisible3-20mm.xlsx"
  frocData <- DfReadDataFile(
    fileName, 
    format = "JAFROC")
  rocData <- DfFroc2Roc(frocData)
  rm(frocData)
}

UtilOutputReport(
  dataset = rocData,
  fom = "Wilcoxon",
  method = "DBMH", 
  reportFormat = "xlsx",
  reportFile = "DBMH.xlsx", 
  showWarnings = FALSE)
UtilOutputReport(
  dataset = rocData,
  fom = "Wilcoxon",
  method = "ORH", 
  reportFormat = "xlsx",
  reportFile = "ORH.xlsx", 
  showWarnings = FALSE)
