rm(list = ls())
library(RJafroc)


################################################################################
# Following lines illustrate the examples on page 16 of the manuscript
################################################################################
str(rocData)
str(frocData)
str(roiData)


################################################################################
# Following lines illustrate the examples on page 18 of the manuscript
################################################################################
rocXlsx <- "http://www.devchakraborty.com/RocData/rocData.xlsx"
rocLrc <- "http://www.devchakraborty.com/RocData/rocData.lrc"
rocCsv <- "http://www.devchakraborty.com/RocData/rocData.csv"
rocImrmc <- "http://www.devchakraborty.com/RocData/rocData.imrmc"
frocXlsx <- "http://www.devchakraborty.com/FrocData/frocData.xlsx"
roiXlsx <- "http://www.devchakraborty.com/RoiData/roiData.xlsx"

fullName <- rocXlsx
download.file(url = fullName, basename(fullName), mode = "wb")
RocDataXlsx<- ReadDataFile(fileName = basename(fullName))

fullName <- rocLrc
download.file(url = fullName, basename(fullName))
RocDataLrc<- ReadDataFile(fileName = basename(fullName), format = "MRMC")

fullName <- rocCsv
download.file(url = fullName, basename(fullName))
RocDataCsv<- ReadDataFile(fileName = basename(fullName), format = "MRMC")

fullName <- rocImrmc
download.file(url = fullName, basename(fullName))
RocDataImrmc<- ReadDataFile(fileName = basename(fullName), format = "iMRMC")

fullName <- frocXlsx
download.file(url = fullName, basename(fullName), mode = "wb")
FrocDataXlsx<- ReadDataFile(fileName = basename(fullName))

fullName <- roiXlsx
download.file(url = fullName, basename(fullName), mode = "wb")
RoiDataXlsx<- ReadDataFile(fileName = basename(fullName))


################################################################################
# Following lines illustrate the examples on page 19 of the manuscript
################################################################################
retDbmRoc  <- DBMHAnalysis(rocData, fom = "Wilcoxon") 
print(retDbmRoc)

retORRoc  <- ORHAnalysis(rocData, fom = "Wilcoxon") 
print(retORRoc)
CovOR <- retORRoc$varComp
cov1 <- CovOR$varCov[3]
cov2 <- CovOR$varCov[4]
cov3 <- CovOR$varCov[5]
varEps <- CovOR$varCov[6]
msTR <- retORRoc$msTR
msT <- retORRoc$msT

print(CovOR)


################################################################################
# Following lines illustrate the examples on page 21 of the manuscript
################################################################################
retDbm  <- DBMHAnalysis(rocData, fom = "Wilcoxon") 
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
varYTR <- retDbm$varComp$varComp[3]
varYTC <- retDbm$varComp$varComp[4]
varYEps <- retDbm$varComp$varComp[6]

for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, varYTR, varYTC, varYEps, 
                          effectSize = effectSize) 
  message("# of readers = ", J, ", estimated # of cases = ", ret$K, "\n",
          "predicted power = ", signif(ret$power, 4), "\n")
}


################################################################################
# Following lines illustrate the examples on page 24 of the manuscript
################################################################################
retOR  <- ORHAnalysis(rocData, fom = "Wilcoxon") 
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
CovOR <- retOR$varComp
cov1 <- CovOR$varCov[3]
cov2 <- CovOR$varCov[4]
cov3 <- CovOR$varCov[5]
varErrOR <- CovOR$varCov[6]
msTR <- retOR$msTR
KStar <- length(rocData$NL[1,1,,1])
for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, cov1 = cov1, cov2 = cov2, cov3 = cov3, 
                          varEps = varErrOR, msTR = msTR, KStar = KStar,
                          effectSize = effectSize) 
  message("# of readers = ", J, ", estimated # of cases = ", ret$K, "\n",
          "predicted power = ", signif(ret$power, 4), "\n")
}


################################################################################
# Following lines illustrate the examples on page 25-28 of the manuscript
# Results are summarized in Tables 7-9
################################################################################
## default JAFROC analysis, wJAFROC FOM is assumed
retDbmwJafroc  <- DBMHAnalysis(frocData)
print(retDbmwJafroc)

retDbmwJafroc1  <- DBMHAnalysis(frocData, fom = "wJAFROC1")
print(retDbmwJafroc1)

retDbmJafroc  <- DBMHAnalysis(frocData, fom = "JAFROC")
print(retDbmJafroc)

## JAFROC1 FOM (use only if there are no non-diseased cases)
retDbmJafroc1  <- DBMHAnalysis(frocData, fom = "JAFROC1")
print(retDbmJafroc1)

# following three examples are for ROC data inferred from FROC data using different methods
retDbmHrAuc  <- DBMHAnalysis(frocData, fom = "HrAuc") # highest rating inferred ROC
retDbmSongA1  <- DBMHAnalysis(frocData, fom = "SongA1")
retDbmSongA2  <- DBMHAnalysis(frocData, fom = "SongA2") # warning: this is very time consuming

# ROI example
retDbmRoi  <- DBMHAnalysis(roiData, fom = "ROI")


################################################################################
# Following lines illustrate the examples on page 28-29 of the manuscript
################################################################################
OutputReport(dataset = rocData, method = "DBMH", fom = "Wilcoxon",
             dataDscrpt = "MyROCData",  showWarnings = FALSE) # showWarnings is set to FALSE as otherwise it will prompt 
                                                              # user "should an existing file be overwritten?"
OutputReport(dataset = rocData, method = "DBMH", fom = "Wilcoxon",
             reportFile = "MyROCDataAnalysis.txt",  showWarnings = FALSE)
OutputReport(dataset = rocData, method = "ORH", fom = "Wilcoxon",
             showWarnings = FALSE)# Following is an example of a wrong usage of the function
OutputReport(dataset = frocData, method = "DBMH", fom = "Wilcoxon",
             showWarnings = FALSE) # ERROR!
OutputReport(dataset = frocData, method = "ORH",
             showWarnings = FALSE) # default fom is wJAFROC
OutputReport(dataset = frocData, method = "DBMH", fom = "HrAuc",
             showWarnings = FALSE)
OutputReport(dataset = roiData, method = "ORH", fom = "ROI",
             showWarnings = FALSE)
# Read data file directly
OutputReport("rocData.xlsx", format = "JAFROC", method = "DBMH",
             fom = "Wilcoxon", dataDscrpt = "MyROC2Data",
             showWarnings = FALSE)


################################################################################
# Following lines illustrate the examples on page 29-30 of the manuscript
# Figure 1(a) and (b)
################################################################################
plotM <- c(1:2)
plotR <- c(1:5)
plotROC <- EmpiricalOpCharac(data = rocData, trts = plotM, 
                             rdrs = plotR, opChType = "ROC")
plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:5),c(1:5))
plotRocAvg <- EmpiricalOpCharac(dataset = rocData, trts = plotMAvg, 
                                rdrs = plotRAvg, opChType = "ROC")
print(plotROC$ROCPlot)
print(plotRocAvg$ROCPlot)


################################################################################
# Following lines illustrate the examples on page 30-31 of the manuscript
# Figure 2(a)-(d)
################################################################################
plotM <- c(1:2)
plotR <- c(1:4)
plotROC <- EmpiricalOpCharac(data = frocData, trts = plotM, 
                             rdrs = plotR, opChType = "ROC")
plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
plotRocAvg <- EmpiricalOpCharac(data = frocData, trts = plotMAvg, 
                                rdrs = plotRAvg, opChType = "ROC")
plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
plotAFROC <- EmpiricalOpCharac(data = frocData, trts = plotMAvg, 
                               rdrs = plotRAvg, opChType = "AFROC")
plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
plotFROC <- EmpiricalOpCharac(data = frocData, trts = plotMAvg, 
                              rdrs = plotRAvg, opChType = "FROC")
print(plotROC$ROCPlot)
print(plotRocAvg$ROCPlot)
print(plotAFROC$AFROCPlot)
print(plotFROC$FROCPlot)

