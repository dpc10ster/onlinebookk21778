install.packages(pkgs = "~/Desktop/RJafroc_0.0.1.tar.gz", repo = NULL, type = "source")
rm(list = ls())
#system("rm *.txt")
library(RJafroc)

##########################################################################################
### read data files from website
rocXlsx <- "http://www.devchakraborty.com/RocData/rocData.xlsx"
rocLrc <- "http://www.devchakraborty.com/RocData/rocData.lrc"
rocCsv <- "http://www.devchakraborty.com/RocData/rocData.csv"
rocImrmc <- "http://www.devchakraborty.com/RocData/rocData.imrmc"
frocXlsx <- "http://www.devchakraborty.com/FrocData/frocData.xlsx"
roiXlsx <- "http://www.devchakraborty.com/RoiData/roiData.xlsx"

fullName <- rocXlsx
download.file(url = fullName,basename(fullName))
RocDataXlsx<- ReadDataFile(fileName = basename(fullName))

fullName <- rocLrc
download.file(url = fullName,basename(fullName))
RocDataLrc<- ReadDataFile(fileName = basename(fullName), format = "MRMC")

fullName <- rocCsv
download.file(url = fullName,basename(fullName))
RocDataCsv<- ReadDataFile(fileName = basename(fullName), format = "MRMC")

fullName <- rocImrmc
download.file(url = fullName,basename(fullName))
RocDataImrmc<- ReadDataFile(fileName = basename(fullName), format = "iMRMC")

fullName <- frocXlsx
download.file(url = fullName,basename(fullName))
FrocDataXlsx<- ReadDataFile(fileName = basename(fullName))

fullName <- roiXlsx
download.file(url = fullName,basename(fullName))
RoiDataXlsx<- ReadDataFile(fileName = basename(fullName))


##########################################################################################
### DBMH analysis
retDbmRoc  <- DBMHAnalysis(rocData, fom = "Wilcoxon") 
varYTR <- retDbmRoc$varComp$varComp[3]
varYTC <- retDbmRoc$varComp$varComp[4]
varYEps <- retDbmRoc$varComp$varComp[6]
effectSize <- retDbmRoc$ciDiffTrtRRRC$Estimate

##########################################################################################
### ORH analysis
retORRoc  <- ORHAnalysis(rocData, fom = "Wilcoxon") 
CovOR <- retORRoc$varComp
cov1 <- CovOR$varCov[3]
cov2 <- CovOR$varCov[4]
cov3 <- CovOR$varCov[5]
varEps <- CovOR$varCov[6]
msTR <- retORRoc$msTR
msT <- retORRoc$msT

##########################################################################################
### sample size
for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, varYTR, varYTC, varYEps, effectSize = abs(effectSize)) 
  cat("# of readers = ", J, "estimated # of cases = ", ret$K, ", predicted power = ", ret$power, "\n")
}

KStar <- length(RocDataXlsx$NL[1,1,,1])
for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, cov1 = cov1, cov2 = cov2, cov3 = cov3, varEps = varEps, msTR = msTR, KStar=KStar, effectSize = abs(effectSize)) 
  cat("# of readers = ", J, "estimated # of cases = ", ret$K, ", predicted power = ", ret$power, "\n")
}

SampleSizeGivenJ(J = 6, cov1 = cov1, cov2 = cov2, cov3 = cov3, varEps = varEps, msTR = msTR, KStar=KStar, effectSize = abs(effectSize))

J <- 6
SampleSizeGivenJ(J, varYTR, varYTC, varYEps, alpha = 0.05,
                 effectSize = abs(effectSize), desiredPower = 0.8, randomOption = "ALL")

K <- 186
PowerGivenJK(J, K, varYTR, varYTC, varYEps, alpha = 0.05,
             effectSize = abs(effectSize), randomOption = "ALL")

PowerTable(data = rocData, alpha = 0.05, effectSize = abs(effectSize), desiredPower = 0.8, randomOption = "ALL")


##########################################################################################
### output reports
rep1 <- OutputReport(data = RocDataXlsx, dataDscrpt = "ROC Data", method = "ORH", fom = "Wilcoxon",
             covEstMethod = "Jackknife", showWarnings = FALSE)
rep2 <- OutputReport(dataset = rocData, method = "DBMH", fom = "Wilcoxon", dataDscrpt = "MyROCData", showWarnings = FALSE)
rep3 <- OutputReport(dataset = rocData, method = "DBMH", fom = "Wilcoxon", reportFile = "MyROCDataAnalysis.txt", showWarnings = FALSE)
##rep4 <- OutputReport(dataset = rocData, method = "ORH", fom = "Wilcoxon", showWarnings = FALSE) ## error
##rep5 <- OutputReport(dataset = frocData, method = "DBMH", fom = "Wilcoxon", showWarnings = FALSE)
rep6 <- OutputReport(dataset = frocData, method = "ORH", showWarnings = FALSE) # default fom is wJAFROC
rep7 <- OutputReport(dataset = frocData, method = "DBMH", fom = "HrAuc", showWarnings = FALSE)
rep8 <- OutputReport(dataset = roiData, method = "ORH", fom = "ROI", showWarnings = FALSE)
rep9 <- OutputReport("rocData.xlsx", format = "JAFROC", method = "DBMH", fom = "Wilcoxon", dataDscrpt = "MyROC2Data", showWarnings = FALSE)


##########################################################################################
### save file in alternate formats
SaveDataFile(dataset = rocData, fileName = "rocData2.xlsx", format = "JAFROC")
SaveDataFile(dataset = rocData, fileName = "rocData2.csv", format = "MRMC")
SaveDataFile(dataset = rocData, fileName = "rocData2.lrc", format = "MRMC", dataDscrpt = "ExampleROCdata")
SaveDataFile(dataset = rocData, fileName = "rocData2.txt", format = "MRMC", dataDscrpt = "ExampleROCdata2")
SaveDataFile(dataset = rocData, fileName = "rocData.imrmc", format = "iMRMC", dataDscrpt = "ExampleROCdata3")


##########################################################################################
### plots
plotM <- c(1:2)
plotR <- c(1:5)
plotROC <- EmpiricalOpCharac(data = RocDataXlsx, modalities = plotM, readers = plotR, OpChType = "ROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:5),c(1:5))
plotRocAvg <- EmpiricalOpCharac(dataset = RocDataXlsx, modalities = plotMAvg, readers = plotRAvg, OpChType = "ROC")
