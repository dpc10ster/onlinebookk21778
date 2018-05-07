rm(list = ls()) #mainSsOrh.R
library(RJafroc)

OptimisticScenario <-  FALSE
fileName <- "VanDyke.lrc" # fileName <- "Franken1.lrc"
cat("File name = ", fileName, "\n")
if (fileName == "VanDyke.lrc") dataset <- dataset02 else dataset <- dataset03

KStar <- length(dataset$NL[1,1,,1])
retOR <- StSignificanceTesting(
  dataset, 
  FOM = "Wilcoxon", 
  method = "ORH") 
effectSize <- retOR$ciDiffTrtRRRC$Estimate
varCompOR <- retOR$varComp
varTR <- varCompOR$varCov[2]
cov1 <- varCompOR$varCov[3]
cov2 <- varCompOR$varCov[4]
cov3 <- varCompOR$varCov[5]
varEps <- varCompOR$varCov[6]
varTR <- max(varTR,0)
effectSize <-retOR$ciDiffTrtRRRC$Estimate
cat("observed effect size =", effectSize, "\n")
sigma <- (retOR$ciDiffTrtRRRC$`CI Upper`-retOR$ciDiffTrtRRRC$`CI Lower`)/4
if (OptimisticScenario == TRUE) {
  if (fileName == "VanDyke.lrc") {
    effectSize <- effectSize -2*sigma
  }
  if (fileName == "Franken1.lrc") {
    effectSize <- effectSize +2*sigma
  }
}
cat("p-value = ", retOR$pRRRC, ", postulated effectSize = ",
    effectSize, ", CI Lower =", retOR$ciDiffTrtRRRC$`CI Lower`,
    ", CI Upper =", retOR$ciDiffTrtRRRC$`CI Upper`, "\n")
powTab <- SsPowerTable(
  effectSize = effectSize, 
  method = "ORH",
  KStar = KStar,
  varTR = varTR, 
  cov1 = cov1, 
  cov2 = cov2, 
  cov3 = cov3, 
  varEps = varEps) 
print(powTab)