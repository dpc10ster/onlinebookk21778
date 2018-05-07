# MainAnalysisNicoMethod.R
rm(list = ls()) 
library(RJafroc)

cat("Repeating CAD ratings 9 times to simulate one modality with 9 identical CAD readers \n9 radiologists form the other modality, followed by ORH analysis")
cat("\nHupse Karssemeijer radiologist only data:\n")
dataNico <- dataset09
NL <- dataNico$NL;LL <- dataNico$LL
K <- length(NL[1,1,,1])
K2 <- length(LL[1,1,,1])
K1 <- K - K2
FP <- NL[1,,1:K1,1]
TP <- LL[1,,,1]

combinedNL <- array(-Inf, dim=c(2,9,K,1))
for (j in 1:9){
  combinedNL[1,j,1:K1,1] <- FP[1,]
}
combinedNL[2,,1:K1,1] <- FP[2:10,]

combinedLL <- array(-Inf, dim=c(2,9,K2,1))
for (j in 1:9){
  combinedLL[1,j,,1] <- TP[1,]
}
combinedLL[2,,,1] <- TP[2:10,]

dataCombined <- dataNico
dataCombined$NL <- combinedNL
dataCombined$LL <- combinedLL
modalityID <- as.character(seq(1,2))
readerID <- as.character(seq(1,9))
dataCombined$modalityID <- modalityID
dataCombined$readerID <- readerID
stats <- StSignificanceTesting(dataCombined,fom = "Wilcoxon", method = "ORH", option = "RRRC")
cat("OR variance components\n")
print(stats$varComp)
cat("F statistic = ", stats$fRRRC, "\n")
cat("ORH ddf (ndf = 1)", stats$ddfRRRC, "\n")
cat("OR p-value = ", stats$pRRRC, "\n")
