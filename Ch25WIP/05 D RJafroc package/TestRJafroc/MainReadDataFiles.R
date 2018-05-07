# the following line needs to be run only once on your machine, to install the software
# install.packages(pkgs = "~/Desktop/RJafroc_0.0.1.tar.gz", repo = NULL, type = "source")
rm(list = ls())
system("rm *.txt");system("rm *.xlsx");system("rm *.lrc");system("rm *.csv");system("rm *.imrmc")
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

