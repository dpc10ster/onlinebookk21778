rm(list = ls()) #mainwAFROCPowerORH.R
library(ggplot2)
library(RJafroc)

stop("fix me")
# included datasets
fileNames <-  c("TONY", "VD", "FR", "FED", "JT", "MAG", "OPT", "PEN", "NICO",
                "RUS", "DOB1", "DOB2", "DOB3", "FZR")
fileName <- fileNames[fileNames == "FED"]
cat("fileName = ", fileName,"\n")
frocData <- loadDataFile(fileName, pathName)
retFileName <- paste0("ANALYZED/RSM-PROPROC-CBM/", "saveRetRoc", fileName)

if (!file.exists((retFileName))){
  I <- length(frocData$modalityID)
  J <- length(frocData$readerID)
  s <- 1
  retSmRoc <- list() # XZ change line 20 to 30
  for (i in 1:I){
    for (j in 1:J){
      cat("i = ", i, ", j = ", j, "\n")
      tempRoc <- FitRsmRoc(frocData, i, j) # fit to RSM
      nLesDistr <- lesionDistribution(frocData)
      retSmRoc[[s]] <- as.list(c(tempRoc, list(nLesDistr = nLesDistr, i = i, j = j)))
      s <- s + 1
    }
  }
  save(retSmRoc, file = retFileName)
} else {
  load(retFileName) # loads object retSmRoc, i.e., ROC data previously analyzed by RSM; this has parameter values
}

i1 <- 1;i2 <- 2 # FED data has 5 modalities; we choose to analyze the first two
cat("NH i1 = ", i1, ", i2 = ", i2, "\n")
selectJ <- c(1, 3, 4, 5) # reader 2 is missing in this dataset
J <- length(frocData$readerID)
KStar <- dim(frocData$NL)[3]

mu1 <- rep(NA, J);mu2 <- rep(NA, J);nu1 <- rep(NA, J);nu2 <- rep(NA, J);lambda1 <- rep(NA, J);lambda2 <- rep(NA, J)

S <- length(retSmRoc)
for (s in 1:S){
  if (retSmRoc[[s]]$mu$Reader %in% selectJ){
    if (retSmRoc[[s]]$mu$Modality == i1){
      mu1[retSmRoc[[s]]$mu$Reader == selectJ] <- as.numeric(retSmRoc[[s]]$mu$mu)
      lambda1[retSmRoc[[s]]$mu$Reader == selectJ] <- as.numeric(retSmRoc[[s]]$lambda$lambda)
      nu1[retSmRoc[[s]]$mu$Reader == selectJ] <- as.numeric(retSmRoc[[s]]$nu$nu)
    }else if (retSmRoc[[s]]$mu$Modality == i2){
      mu2[retSmRoc[[s]]$mu$Reader == selectJ] <- as.numeric(retSmRoc[[s]]$mu$mu)
      lambda2[retSmRoc[[s]]$mu$Reader == selectJ] <- as.numeric(retSmRoc[[s]]$lambda$lambda)
      nu2[retSmRoc[[s]]$mu$Reader == selectJ] <- as.numeric(retSmRoc[[s]]$nu$nu)
    }
  }
}

mu <- rbind(mu1, mu2);lambda <- rbind(lambda1, lambda2);nu <- rbind(nu1, nu2)
muMed <- median(mu) # instead of average, use median to get representative value over whole dataset
nuMed <- median(nu) # do:
lambdaMed <- median(lambda) # do:

# construct lesion weights, assuming equally weighted lesions
nLesDistr <- retSmRoc[[1]]$nLesDistr
lesionWeights <- matrix(-Inf, nrow = nrow(nLesDistr), ncol = nrow(nLesDistr))
for (l in 1:nrow(nLesDistr)){
  nLes <- nLesDistr[l, 1]
  lesionWeights[l, 1:nLes] <- 1/nLes
}
# calculate NH values for ROC-AUC and wAFROC-AUC
aucRocNH <- PlotRsmOperatingCharacteristics(muMed, lambdaMed, nuMed, 
                                        lesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "ROC")$aucROC
aucAfrocNH <- PlotRsmOperatingCharacteristics(muMed, lambdaMed, nuMed, 
                                          lesionDistribution = nLesDistr, lesionWeights = lesionWeights, type = "wAFROC")$aucwAFROC

# following code calculates ROC-ES and wAFROC-ES
deltaMu <- seq(0.01, 0.2, 0.01) # values of deltaMu to scan below
esRoc <- array(dim = length(deltaMu));eswAfroc <- array(dim = length(deltaMu))
for (i in 1:length(deltaMu)) {
  esRoc[i] <- PlotRsmOperatingCharacteristics(muMed + deltaMu[i], lambdaMed, nuMed, lesionDistribution = 
                                            nLesDistr, lesionWeights = lesionWeights, type = "ROC")$aucROC - aucRocNH
  eswAfroc[i] <- PlotRsmOperatingCharacteristics(muMed+ deltaMu[i], lambdaMed, nuMed, lesionDistribution = 
                                               nLesDistr, lesionWeights = lesionWeights, type = "wAFROC")$aucwAFROC - aucAfrocNH
  #cat("ES ROC, wAFROC = ", esRoc[i], eswAfroc[i],"\n")
}
#cat("\n")

a<-lm(eswAfroc~-1+esRoc) # fit values to straight line thru origin
effectSizeROC <- seq(0.01, 0.1, 0.01)
effectSizewAFROC <- effectSizeROC*a$coefficients[1]

JTest <- 5;KTest <- 100
varCompROC <- StSignificanceTesting(frocData, fom = "HrAuc", method = "ORH", option = "RRRC")$varComp
varCompwAFROC <- StSignificanceTesting(frocData, fom = "wAFROC", method = "ORH", option = "RRRC")$varComp

cat("JTest = ", JTest, "KTest = ", KTest, "\n")
powerROC <- array(dim = length(effectSizeROC));powerwAFROC <- array(dim = length(effectSizeROC))
for (i in 1:length(effectSizeROC)) {
  cov1 <- varCompROC$varCov[3]
  cov2 <- varCompROC$varCov[4]
  cov3 <- varCompROC$varCov[5]
  varTR <- varCompROC$varCov[2]
  varEps <- varCompROC$varCov[6]
  powerROC[i] <- SsPowerGivenJK(JTest, KTest, alpha = 0.05, effectSize = effectSizeROC[i], option = "RRRC",
                                method = "ORH", cov1 = cov1, cov2 = cov2, cov3 = cov3, varTR = varTR, 
                                varEps = varEps, KStar = KStar)$powerRRRC
  
  cov1 <- varCompwAFROC$varCov[3]
  cov2 <- varCompwAFROC$varCov[4]
  cov3 <- varCompwAFROC$varCov[5]
  varTR <- varCompwAFROC$varCov[2]
  varEps <- varCompwAFROC$varCov[6]
  powerwAFROC[i] <- SsPowerGivenJK(JTest, KTest, alpha = 0.05, effectSize = effectSizewAFROC[i], option = "RRRC",
                                   method = "ORH", cov1 = cov1, cov2 = cov2, cov3 = cov3, varTR = varTR, 
                                   varEps = varEps, KStar = KStar)$powerRRRC
  cat("ROC effect size = ,", effectSizeROC[i], "wAFROC effect size = ,", effectSizewAFROC[i], 
      ", Statistical power ROC, wAFROC:", powerROC[i], ",", powerwAFROC[i], "\n")
}

df <- data.frame(esRoc = esRoc, eswAfroc = eswAfroc)
p <- ggplot(data = df, aes(x = esRoc, y = eswAfroc)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_point()
print(p)

df <- data.frame(powerROC = powerROC, powerwAFROC = powerwAFROC)
p <- ggplot(mapping = aes(x = powerROC, y = powerwAFROC)) +
  geom_line(data = df, size = 0.5)
print(p)
