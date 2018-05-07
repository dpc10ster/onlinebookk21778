###################################################################
#ReadDataFile
ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/RocData/rocData.xlsx"
download.file(url = fullName,basename(fullName))
RocDataXlsx<- ReadDataFile(basename(fullName))
proc.time() - ptm
#user  system elapsed 
#2.074   0.079   1.307 

ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/RocData/rocData.csv"
download.file(url = fullName,basename(fullName))
RocDataXlsx<- ReadDataFile(basename(fullName))
proc.time() - ptm
# user  system elapsed 
# 0.011   0.010   0.153 

ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/RocData/rocData.imrmc"
download.file(url = fullName,basename(fullName))
RocDataXlsx<- ReadDataFile(basename(fullName))
proc.time() - ptm
# user  system elapsed 
# 0.010   0.003   0.200

ptm <- proc.time()
fullName <- "http://www.devchakraborty.com/FrocData/frocData.xlsx"
download.file(url = fullName,basename(fullName))
temp <- ReadDataFile(basename(fullName))
proc.time() - ptm
# user  system elapsed 
# 2.259   0.098   1.325 

ptm <- proc.time()
roiXlsx <- "http://www.devchakraborty.com/RoiData/roiData.xlsx"
fullName <- roiXlsx
download.file(url = fullName, basename(fullName))
temp <- ReadDataFile(basename(fullName))
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
retDbmROI  <- DBMHAnalysis(roiData, fom = "ROI")
print(retDbmROI)
proc.time() - ptm
# user  system elapsed 
# 2.103   0.046   2.147 

ptm <- proc.time()
retDbmwJafroc1  <- DBMHAnalysis(frocData, fom = "wJAFROC1")
print(retDbmwJafroc1)
proc.time() - ptm
# user  system elapsed 
# 4.010   0.135   4.144 

ptm <- proc.time()
retDbmwJAFROC  <- DBMHAnalysis(frocData) 
print(retDbmwJAFROC)
proc.time() - ptm
# user  system elapsed 
# 1.489   0.086   1.109



###################################################################
# reports
ptm <- proc.time()
OutputReport(data = rocData, dataDscrpt = "ROC Data", method = "ORH", fom = "Wilcoxon",
             covEstMethod = "Jackknife")
proc.time() - ptm


ptm <- proc.time()
OutputReport("rocData.xlsx", format = "JAFROC", method = "DBMH", fom = "Wilcoxon")
proc.time() - ptm




