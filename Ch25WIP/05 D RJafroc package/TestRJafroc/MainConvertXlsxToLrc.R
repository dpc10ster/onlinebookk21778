rm(list = ls()) # delete all existing variables
#system("rm *.txt");system("rm *.xlsx");system("rm *.lrc");system("rm *.csv");system("rm *.imrmc")
library(RJafroc)
source("ConvertFROC2HRatings.R")

## ret <- ReadDataFile(fileName = "./DataSets/FedericaAll.xlsx")
## ret <- ReadDataFile(fileName = "./DataSets/PET-CT (10AR).xlsx")
ret <- ReadDataFile(fileName = "./DataSets/VanDykeData.xlsx")
NL <- ret$NL;LL <- ret$LL
K2 <- length(LL[1,1,,1])
K1 <- length(NL[1,1,,1])-K2
I <- length(LL[,1,1,1])
J <- length(LL[1,,1,1])

ret2 <- ConvertFROC2HRatings (NL, LL )
FP <- array(dim=c(I,J,K1+K2,1))
FP[,,1:K1,1] <- ret2$Hratings[,,1:K1]
ret$NL <- FP 

TP <- array(dim=c(I,J,K2,1))
TP[,,,1] <- ret2$Hratings[,,(K1+1):(K1+K2)]
ret$LL <- TP 

ret$lesionNum <- rep(1,K2)
ret$dataType <- "ROC"

lesionID <- array(1, dim = c(K2,1))
ret$lesionID <- lesionID

lesionWeight <- array(1, dim = c(K2,1))
ret$lesionWeight <- lesionWeight

## SaveDataFile(dataset = ret,fileName = "./DataSets/FedericaAll.lrc", format = "MRMC")
## SaveDataFile(dataset = ret,fileName = "./DataSets/PET-CT (10AR).lrc", format = "MRMC")
SaveDataFile(dataset = ret,fileName = "./VanDykeData.lrc", format = "MRMC")

