rm(list = ls()) # mainDataSummary.R # freeze lines
library(RJafroc)
source("loadDataFile.R")

pathName <- "../../../06 E Online Appendices/E24 Datasets/"
# included datasets
fileNames <-  c("TONY", "VD", "FR", "FED", "JT", "MAG", "OPT", "PEN", "NICO",
                "RUS", "DOB1", "DOB2", "DOB3", "FZR")
descriptions <- c(
  "Digital breast tomosynthesis vs. mammography",
  "Cine vs. SE MRI for aortic dissection",
  "Digital vs. analog pediatric chest",
  "Image processing in mammography: FROC",
  "Nodule detection in an thorax CT phantom",
  "Tomosynthesis Vs. Radiography Pulmonary Nodules",
  "Calcification detection in digital mammography",
  "Image compression in mammography",
  "Standalone CAD vs. radiologists mammography",
  "Lesion detection in digital mammography",
  "Tomosynthesis, Dual-Energy & Conventional Chest", 
  "Tomosynthesis, Dual-Energy & Conventional Chest",
  "Tomosynthesis, Dual-Energy & Conventional Chest",
  "Image processing in mammography: ROC"
)

fileName <- fileNames[fileNames == "TONY"]
for (f in 1:length(fileNames)){
  frocData <- loadDataFile(fileNames[f], pathName)
  rocData <- DfFroc2HrRoc(frocData)
  I <- length(rocData$modalityID)
  J <- length(rocData$readerID)
  K <- dim(rocData$NL)[3]
  K2 <- dim(rocData$LL)[3]
  K1 <- K - K2
  cat(fileNames[f], I, J, K1, K2, K, descriptions[f], "\n")
}
