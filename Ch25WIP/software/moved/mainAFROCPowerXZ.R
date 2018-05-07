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

mu1 <- mean(mu1)
mu2 <- mean(mu2)

nLesDistr <- retSmRoc[[1]]$nLesDistr
lesionWeights <- matrix(-Inf, nrow = nrow(nLesDistr), ncol = nrow(nLesDistr))
for (l in 1:nrow(nLesDistr)){
  nLes <- nLesDistr[l, 1]
  lesionWeights[l, 1:nLes] <- 1/nLes
}

aucROC1 <- rep(NA, J)
aucROC2 <- aucROC1
for (j in 1:J){
  aucROC1[j] <- OperatingCharacteristics(mu1, lambda1[j], nu1[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "ROC")$aucROC
  aucROC2[j] <- OperatingCharacteristics(mu2, lambda2[j], nu2[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "ROC")$aucROC
}

effectSizeROC <- mean(aucROC2) - mean(aucROC1)
deltaMu <- mu2 * 0.1
while(effectSizeROC < 0.05){
  mu2 <- mu2 + deltaMu
  for (j in 1:J){
    aucROC2[j] <- OperatingCharacteristics(mu2, lambda2[j], nu2[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "ROC")$aucROC
  }
  effectSizeROC <- mean(aucROC2) - mean(aucROC1)
}


aucAFROC1 <- rep(NA, J)
aucAFROC2 <- aucAFROC1
for (j in 1:J){
  aucAFROC1[j] <- OperatingCharacteristics(mu1, lambda1[j], nu1[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "AFROC")$aucAFROC
  aucAFROC2[j] <- OperatingCharacteristics(mu2, lambda2[j], nu2[j], pmfLesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "AFROC")$aucAFROC
}
effectSizeAFROC <- mean(aucAFROC2) - mean(aucAFROC1)

cat("aucROC1: ", aucROC1, "\n")
cat("aucAFROC1: ", aucAFROC1, "\n")

cat("aucROC2: ", aucROC2, "\n")
cat("aucAFROC2: ", aucAFROC2, "\n")

cat("mean(aucROC1):", mean(aucROC1),"\n")
cat("mean(aucAFROC1):", mean(aucAFROC1),"\n")

cat("mean(aucROC2):", mean(aucROC2),"\n")
cat("mean(aucAFROC2):", mean(aucAFROC2),"\n")

stop("temp....")
JTest <- 5
KTest <- 135
varCompROC <- DBMHAnalysis(frocData, fom = "HrAuc", option = "RRRC")$varComp
varYTR <- varCompROC$varComp[3]
varYTC <- varCompROC$varComp[4]
varYEps <- varCompROC$varComp[6]
powerROC <- PowerGivenJK(JTest, KTest, varYTR = varYTR, varYTC = varYTC, varYEps = varYEps, effectSize = effectSizeROC)
cat("Statistical power of ROC:", powerROC, "\n")

varCompAFROC <- DBMHAnalysis(frocData, fom = "AFROC", option = "RRRC")$varComp
varYTR <- varCompAFROC$varComp[3]
varYTC <- varCompAFROC$varComp[4]
varYEps <- varCompAFROC$varComp[6]
effectSizeAFROC <- effectSizeAFROC / 2
powerAFROC <- PowerGivenJK(JTest, KTest, varYTR = varYTR, varYTC = varYTC, varYEps = varYEps, effectSize = effectSizeAFROC)
cat("Statistical power of AFROC:", powerAFROC, "\n")

