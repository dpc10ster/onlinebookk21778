rm(list = ls()) # delete all existing variables
system("rm *.txt");system("rm *.xlsx");system("rm *.lrc");system("rm *.csv");system("rm *.imrmc")
install.packages(pkgs = "~/Desktop/RJafroc_0.0.2.tar.gz", repo = NULL, type = "source")
library(RJafroc)

RJafrocGui()

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
effectSize <- retORRoc$ciDiffTrtRRRC$Estimate
CovOR <- retORRoc$varComp
cov1 <- CovOR$varCov[3]
cov2 <- CovOR$varCov[4]
cov3 <- CovOR$varCov[5]
varEps <- CovOR$varCov[6]
msTR <- retORRoc$msTR
msT <- retORRoc$msT

##########################################################################################
### sample size
### sample size using DBMH variance components
for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, varYTR, varYTC, varYEps, effectSize = effectSize) 
  cat("# of rdrs = ", J, "estimated # of cases = ", ret$K, ", predicted power = ", ret$power, "\n")
}

### sample size using ORH variance components
KStar <- length(rocData$NL[1,1,,1])
for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, cov1 = cov1, cov2 = cov2, cov3 = cov3, varEps = varEps, msTR = msTR, KStar=KStar, effectSize = effectSize) 
  cat("# of rdrs = ", J, "estimated # of cases = ", ret$K, ", predicted power = ", ret$power, "\n")
}

K <- 251
PowerGivenJK(6, K, varYTR, varYTC, varYEps, alpha = 0.05,
             effectSize = effectSize, randomOption = "ALL")

PowerTable(data = rocData, alpha = 0.05, effectSize = effectSize, desiredPower = 0.8, randomOption = "READERS")

##########################################################################################
### output reports
rep1 <- OutputReport(data = rocData, dataDscrpt = "ROC Data", method = "ORH", fom = "Wilcoxon",
             covEstMethod = "Jackknife", showWarnings = FALSE)
rep2 <- OutputReport(dataset = rocData, method = "DBMH", fom = "Wilcoxon", dataDscrpt = "MyROCData", showWarnings = FALSE)
rep3 <- OutputReport(dataset = rocData, method = "DBMH", fom = "Wilcoxon", reportFile = "MyROCDataAnalysis.txt", showWarnings = FALSE)
## rep4 <- OutputReport(dataset = frocData, method = "ORH", fom = "SongA2", showWarnings = FALSE) ## very long computation
## rep5 <- OutputReport(dataset = frocData, method = "DBMH", fom = "Wilcoxon", showWarnings = FALSE) ## error
rep6 <- OutputReport(dataset = frocData, method = "ORH", showWarnings = FALSE) # default fom is wJAFROC
rep7 <- OutputReport(dataset = frocData, method = "DBMH", fom = "HrAuc", showWarnings = FALSE)
rep8 <- OutputReport(dataset = roiData, method = "ORH", fom = "ROI", showWarnings = FALSE)
rep9 <- OutputReport(dataset = rocData, format = "JAFROC", method = "DBMH", fom = "Wilcoxon", dataDscrpt = "MyROC2Data", showWarnings = FALSE)

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
plotRoc <- EmpiricalOpCharac(data = rocData, trts = plotM, rdrs = plotR, opChType = "ROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:5),c(1:5))
plotRocAvg <- EmpiricalOpCharac(dataset = rocData, trts = plotMAvg, rdrs = plotRAvg, opChType = "ROC")
