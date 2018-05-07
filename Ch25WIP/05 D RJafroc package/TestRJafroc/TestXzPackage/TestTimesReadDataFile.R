install.packages(pkgs = "~/Desktop/RJafroc_0.0.1.tar.gz", repo = NULL, type = "source")
#install.packages(pkgs = "~/Desktop/RJafroc_0.0-1.tar.gz", repo = NULL, type = "source")

rm(list = ls())
library(RJafroc)

###################################################################
# ReadDataFiles
ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/RocData/rocData.xlsx"
download.file(url = fullName,basename(fullName))
RocDataXlsx<- ReadDataFile(basename(fullName))
proc.time() - ptm
#user  system elapsed 
# 0.402   0.042   0.489  

ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/RocData/rocData.csv"
download.file(url = fullName,basename(fullName))
RocDataCsv<- ReadDataFile(basename(fullName))
proc.time() - ptm
# user  system elapsed 
# 0.011   0.010   0.153 

ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/RocData/rocData.imrmc"
download.file(url = fullName,basename(fullName))
RocDataImrmc<- ReadDataFile(basename(fullName))
proc.time() - ptm
# user  system elapsed 
# 0.010   0.003   0.200

ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/FrocData/frocData.xlsx"
download.file(url = fullName,basename(fullName))
FrocDataXlsx <- ReadDataFile(basename(fullName))
proc.time() - ptm
# user  system elapsed 
# 2.259   0.098   1.325 

ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/RoiData/roiData.xlsx"
download.file(url = fullName, basename(fullName))
RoiDataXlsx <- ReadDataFile(basename(fullName))
proc.time() - ptm
# user  system elapsed 
# 1.489   0.086   1.109 

###################################################################
#analysis
ptm <- proc.time()
retDbmRoc  <- DBMHAnalysis(rocData, fom = "Wilcoxon") 
print(retDbmRoc)
proc.time() - ptm
# user  system elapsed 
# 0.698   0.035   0.730 

ptm <- proc.time()
retORRoc  <- ORHAnalysis(rocData, fom = "Wilcoxon") 
print(retORRoc)
proc.time() - ptm
# user  system elapsed 
# 1.322   0.084   1.406 

ptm <- proc.time()
retDbmSongA1  <- DBMHAnalysis(frocData, fom = "SongA1") 
print(retDbmSongA1)
proc.time() - ptm
# user  system elapsed 
# 22.531   0.047  22.570 

ptm <- proc.time()
retDbmHrAuc  <- DBMHAnalysis(frocData, fom = "HrAuc") 
print(retDbmHrAuc)
proc.time() - ptm
# user  system elapsed 
# 2.404   0.135   2.541 

ptm <- proc.time()
retDbmSongA2  <- DBMHAnalysis(frocData, fom = "SongA2") 
print(retDbmSongA2)
proc.time() - ptm
# user  system elapsed 
# 142.707   0.313 143.006 

ptm <- proc.time()
retDbmwJafroc1  <- DBMHAnalysis(frocData, fom = "wJAFROC1")
print(retDbmwJafroc1)
proc.time() - ptm
# user  system elapsed 
# 4.010   0.135   4.144 

ptm <- proc.time()
retDbmwJAFROC  <- DBMHAnalysis(frocData) # default is weighted JAFROC
print(retDbmwJAFROC)
proc.time() - ptm
# user  system elapsed 
# 1.489   0.086   1.109

ptm <- proc.time()
retDbmJafroc1  <- DBMHAnalysis(frocData, fom = "JAFROC1")
print(retDbmJafroc1)
proc.time() - ptm
# user  system elapsed 
# 3.200   0.230   3.428 

ptm <- proc.time()
retDbmJAFROC  <- DBMHAnalysis(frocData, fom = "JAFROC") 
print(retDbmJAFROC)
proc.time() - ptm
# user  system elapsed 
# 1.884   0.124   2.009 

system.time(retDbmROI  <- DBMHAnalysis(roiData, fom = "ROI"))
print(retDbmROI)
# user  system elapsed 
# 2.103   0.046   2.147 

###################################################################
# Sample size estimates
## Following is an example of sample size calculation with DBM variance components.
ptm <- proc.time()
retDbm <- DBMHAnalysis(data = rocData, fom = "Wilcoxon")
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
varCompDBM <- retDbm$varComp
varYTR <- varCompDBM$varComp[3]
varYTC <- varCompDBM$varComp[4]
varYEps <- varCompDBM$varComp[6]
SampleSizeGivenJ(J = 6, varYTR = varYTR, varYTC = varYTC, varYEps = varYEps, effectSize = effectSize)
proc.time() - ptm
# user  system elapsed 
# 0.670   0.034   0.748 iMac

ptm <- proc.time()
retOR <- ORHAnalysis(data = rocData, fom = "Wilcoxon", covEstMethod = "Jackknife")
effectSize <- retOR$ciDiffTrtRRRC$Estimate
varCompOR <- retOR$varComp
cov1 <- varCompOR$varCov[3]
cov2 <- varCompOR$varCov[4]
cov3 <- varCompOR$varCov[5]
varEps <- varCompOR$varCov[6]
msTR <- retOR$msTR
KStar <- 114
SampleSizeGivenJ(J = 6, cov1 = cov1, cov2 = cov2, cov3 = cov3, varEps= varEps,
                 msTR = msTR, KStar = KStar, effectSize = effectSize)
proc.time() - ptm
# user  system elapsed 
# 1.295   0.090   1.435 iMac

ptm <- proc.time()
retDbm <- DBMHAnalysis(data = rocData, fom = "Wilcoxon")
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
varYTR <- retDbm$varComp$varComp[3]
varYTC <- retDbm$varComp$varComp[4]
varYEps <- retDbm$varComp$varComp[6]
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
for (J in 6:10) {
  ret <- SampleSizeGivenJ(J, varYTR, varYTC, varYEps, effectSize = effectSize)
  message("# of readers = ", J, "estimated # of cases = ", ret$K, ", predicted power = ",
          signif(ret$power,3), "\n")
}
proc.time() - ptm
# user  system elapsed 
# 0.716   0.036   0.751 


# ptm <- proc.time()
# retOR <- ORHAnalysis(data = rocData, fom = "Wilcoxon", covEstMethod = "Bootstrap")
# effectSize <- retOR$ciDiffTrtRRRC$Estimate
# varCompOR <- retOR$varComp
# cov1 <- varCompOR$varCov[3]
# cov2 <- varCompOR$varCov[4]
# cov3 <- varCompOR$varCov[5]
# varEps <- varCompOR$varCov[6]
# msTR <- retOR$msTR
# KStar <- length(rocData$NL[1,1,,1])
# SampleSizeGivenJ(J = 6, cov1 = cov1, cov2 = cov2, cov3 = cov3, varEps= varEps,
#                  msTR = msTR, KStar = KStar, effectSize = effectSize)
# proc.time() - ptm
# # user  system elapsed 
# # 2.338   0.145   2.484 

ptm <- proc.time()
retDbm <- DBMHAnalysis(data = rocData, fom = "Wilcoxon")
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
PowerTable(dataset = rocData, effectSize = effectSize)
proc.time() - ptm
# user  system elapsed 
# 1.632   0.058   1.691 


## Following is an example of sample size calculation with DBM variance componements.
retDbm <- DBMHAnalysis(data = rocData, fom = "Wilcoxon")
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
varCompDBM <- retDbm$varComp
varYTR <- varCompDBM$varComp[3]
varYTC <- varCompDBM$varComp[4]
varYEps <- varCompDBM$varComp[6]
PowerGivenJK(numReaders = 6, numCases = 186, varYTR = varYTR, varYTC = varYTC,
             varYEps = varYEps, effectSize = effectSize)

## Following is an example of sample size calculation with OR variance componements.
retOR <- ORHAnalysis(data = rocData, fom = "Wilcoxon", covEstMethod = "Jackknife")
effectSize <- retOR$ciDiffTrtRRRC$Estimate
varCompOR <- retOR$varComp
cov1 <- varCompOR$varCov[3]
cov2 <- varCompOR$varCov[4]
cov3 <- varCompOR$varCov[5]
varEps <- varCompOR$varCov[6]
KStar <- length(rocData$NL[1,1,,1])
msTR <- retOR$msTR
PowerGivenJK(numReaders = 6, numCases = 186, cov1 = cov1, cov2 = cov2, cov3 = cov3,
             varEps = varEps, KStar = KStar, msTR = msTR, effectSize = effectSize)




###################################################################
# reports
system.time(OutputReport(dataset = rocData, method = "DBMH", fom = "Wilcoxon", dataDscrpt = "MyROCData", showWarnings = FALSE))
# 0.680   0.037   0.718 

system.time(OutputReport(dataset = rocData, method = "DBMH", fom = "Wilcoxon", reportFile = "MyROCDataAnalysis.txt", showWarnings = FALSE)) 
# 0.676   0.039   0.716 

system.time(OutputReport(dataset = rocData, method = "ORH", fom = "Wilcoxon", showWarnings = FALSE))
# 1.271   0.086   1.357 

system.time(OutputReport(dataset = frocData, method = "ORH", showWarnings = FALSE)) # default fom is wJAFROC
# 5.662   0.139   5.797 

system.time(OutputReport(dataset = frocData, method = "DBMH", fom = "HrAuc", showWarnings = FALSE))
# 2.414   0.127   2.540 

system.time(OutputReport(dataset = roiData, method = "ORH", fom = "ROI", showWarnings = FALSE))
# 5.558   0.128   5.683

