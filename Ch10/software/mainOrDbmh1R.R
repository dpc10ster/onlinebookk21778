rm(list = ls()) #mainOrDbmh1R.R
library(RJafroc)
source("Wilcoxon.R")
source("VarCov1Jk.R")
fileName <- "CXRinvisible3-20mm.xlsx"
frocData <- DfReadDataFile(fileName, format = "JAFROC", newExcelFileFormat = FALSE)
rocData <- DfFroc2Roc(frocData)
jSelect <- 1
rocData1R <- DfExtractDataset(rocData, rdrs = jSelect)
cat("data file = ", fileName, "\n")
cat("selected reader = ", jSelect, "\n")
seed <- 1; set.seed(seed)

ret1 <- StSignificanceTesting(
  rocData1R,
  FOM = "Wilcoxon", 
  method = "DBMH", 
  option = "FRRC")
cat("DBMH: F-stat = ", 
    ret1$fFRRC, 
    "\nddf = ", ret1$ddfFRRC, 
    "\nP-val = ", ret1$pFRRC,"\n")

ret2 <- StSignificanceTesting(
  rocData1R,
  FOM = "Wilcoxon", 
  method = "ORH", 
  option = "FRRC")
cat("ORH (Jackknife):  F-stat = ", 
    ret2$fFRRC, "\nddf = ", 
    ret2$ddfFRRC, "\nP-val = ", 
    ret2$pFRRC,"\n")

ret3 <- StSignificanceTesting(
  rocData1R,
  FOM = "Wilcoxon", 
  method = "ORH", 
  option = "FRRC", 
  covEstMethod = "DeLong")
cat("ORH (DeLong):  F-stat = ", 
    ret3$fFRRC, 
    "\nddf = ", ret3$ddfFRRC, "\nP-val = ", 
    ret3$pFRRC,"\n")

ret4 <- StSignificanceTesting(
  rocData1R,
  FOM = "Wilcoxon", 
  method = "ORH", 
  option = "FRRC", 
  covEstMethod = "Bootstrap")
cat("ORH (Bootstrap):  F-stat = ", 
    ret4$fFRRC, "\nddf = ", 
    ret4$ddfFRRC, "\nP-val = ", 
    ret4$pFRRC,"\n")