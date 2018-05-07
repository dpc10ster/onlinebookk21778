install.packages(pkgs = "~/Desktop/RJafroc_0.0-1.tar.gz", repo = NULL, type = "source")
rm(list = ls())
library(RJafroc)

numLesDistr <- rbind(c(1, 35), c(2, 25), c(3, 10))
OperatingCharacteristics(mu = c(2, 3), lambda = c(1, 1.5), nu = c(0.6, 0.8),
                         numLesDistr = numLesDistr, legendPosition = "bottom")

roiData <- ReadDataFile("roiData.xlsx")

retDbmRoi  <- DBMHAnalysis(roiData, fom = "ROI") 
print(retDbmRoi)

OutputReport("roiData.xlsx", format = "JAFROC", method = "DBMH", fom = "ROI")

