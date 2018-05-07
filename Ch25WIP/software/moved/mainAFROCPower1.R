rm(list = ls())

library(RJafroc)
library(ggplot2)
source("OperatingCharacteristics.R")

fileName <- "FED" 
# fileName <- "TONY"
# fileName <- "VD"
# fileName <- "TONY"
# fileName <- "FR"

if (fileName == "FED"){
  frocData <- ReadDataFile("./Datasets/FZ/FedericaAll.xlsx")  
} else if (fileName == "TONY") {
  frocData <- ReadDataFile("./Datasets/ALMLCWBZSZ_Finale_20100402.xlsx")  
} else if (fileName == "VD") {
  frocData <- ReadDataFile("./Datasets/VanDykeData.xlsx")  
} else if (fileName == "JT") {
  frocData <- ReadDataFile("./Datasets/JohnThompsonDatafile.xlsx")  
} else if (fileName == "FR") {
  frocData <- ReadDataFile("./Datasets/franken1.xlsx")  
}
retFileName <- paste0("saveRetRoc", fileName)
load(retFileName)

i1 <- 4
i2 <- 5
selectJ <- c(1, 2, 3, 4)
frocData <- ExtractDataset(frocData, trts = c(i1, i2), rdrs = selectJ)
J <- length(frocData$readerID)
K <- dim(frocData$NL)[3]

mu1 <- rep(NA, J)
mu2 <- mu1
nu1 <- mu1
nu2 <- mu1
lambda1 <- mu1
lambda2 <- mu1

S <- length(retSmRoc)
for (s in 1:S){
  if (retSmRoc[[s]]$j %in% selectJ){
    if (retSmRoc[[s]]$i == i1){
      mu1[retSmRoc[[s]]$j == selectJ] <- as.numeric(retSmRoc[[s]]$mu$mu)
      nu1[retSmRoc[[s]]$j == selectJ] <- as.numeric(retSmRoc[[s]]$nu$nu)
      lambda1[retSmRoc[[s]]$j == selectJ] <- as.numeric(retSmRoc[[s]]$lambda$lambda)
    }else if (retSmRoc[[s]]$i == i2){
      mu2[retSmRoc[[s]]$j == selectJ] <- as.numeric(retSmRoc[[s]]$mu$mu)
      nu2[retSmRoc[[s]]$j == selectJ] <- as.numeric(retSmRoc[[s]]$nu$nu)
      lambda2[retSmRoc[[s]]$j == selectJ] <- as.numeric(retSmRoc[[s]]$lambda$lambda)
    }
  }
}

nLesDistr <- retSmRoc[[1]]$nLesDistr
lesionWeights <- matrix(-Inf, nrow = nrow(nLesDistr), ncol = nrow(nLesDistr))
for (l in 1:nrow(nLesDistr)){
  nLes <- nLesDistr[l, 1]
  lesionWeights[l, 1:nLes] <- 1/nLes
}

Roc1 <- rep(NA, J)
Roc2 <- Roc1
for (j in 1:J){
  Roc1[j] <- OperatingCharacteristics(mu1[j], lambda1[j], nu1[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "ROC")$aucROC
  Roc2[j] <- OperatingCharacteristics(mu2[j], lambda2[j], nu2[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "ROC")$aucROC
}
RocEs <- mean(Roc2) - mean(Roc1)

Afroc1 <- rep(NA, J)
Afroc2 <- Afroc1
for (j in 1:J){
  Afroc1[j] <- OperatingCharacteristics(mu1[j], lambda1[j], nu1[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "AFROC")$aucAFROC
  Afroc2[j] <- OperatingCharacteristics(mu2[j], lambda2[j], nu2[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "AFROC")$aucAFROC
}
AfrocEs <- mean(Afroc2) - mean(Afroc1)

empRoc <- FigureOfMerit(frocData, fom = "HrAuc")
empAfroc <- FigureOfMerit(frocData)
empRocEs <- mean(empRoc[1,]) - mean(empRoc[2,])
empAfrocES <- mean(empAfroc[1,]) - mean(empAfroc[2,])

cat("empRoc1: ", empRoc[1,], "\n")
cat("empAfroc1: ", empAfroc[1,], "\n")
cat("Roc1: ", Roc1, "\n")
cat("Afroc1: ", Afroc1, "\n")
cat("\n")
cat("empRoc2: ", empRoc[2,], "\n")
cat("empAfroc2: ", empAfroc[2,], "\n")
cat("Roc2: ", Roc2, "\n")
cat("Afroc2: ", Afroc2, "\n")

cat("\n")
cat("mean(empRoc[1,]):", mean(empRoc[1,]),"\n")
cat("mean(empAfroc[1,]):", mean(empAfroc[1,]),"\n")
cat("mean(Roc1):", mean(Roc1),"\n")
cat("mean(Afroc1):", mean(Afroc1),"\n")

cat("\n")
cat("mean(empRoc[2,]):", mean(empRoc[2,]),"\n")
cat("mean(empAfroc[2,]):", mean(empAfroc[2,]),"\n")
cat("mean(Roc2):", mean(Roc2),"\n")
cat("mean(Afroc2):", mean(Afroc2),"\n")

cat("\n")
cat("empRocEs:", empRocEs,"\n")
cat("empAfrocES:", empAfrocES,"\n")
cat("RocEs:", RocEs,"\n")
cat("AfrocEs:", AfrocEs,"\n")
