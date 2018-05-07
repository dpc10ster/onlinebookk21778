rm(list = ls()) #mainDBMHBrief.R
library(RJafroc)
library(ggplot2)
#fileName <- "Franken1.lrc"
fileName <- "VanDyke.lrc"
UtilOutputReport(
  fileName, format = "MRMC",
  method = "DBMH", FOM = "Wilcoxon", 
  showWarnings ="FALSE",
  reportFile = "VanDykeOutput.txt")
UtilOutputReport(
  fileName, format = "MRMC", reportFormat = "xlsx",
  method = "DBMH", FOM = "Wilcoxon", 
  showWarnings = "FALSE",
  reportFile = "VanDykeOutput.xlsx")
rocData <- DfReadDataFile(fileName, format = "MRMC")
plot14 <- PlotEmpiricalOperatingCharacteristics(
  rocData, trts = 1, rdrs = 4)
print(plot14$Plot);print(plot14$Points)
plot24 <- PlotEmpiricalOperatingCharacteristics(
  rocData, trts = 2, rdrs = 4)
print(plot24$Plot);print(plot24$Points)
