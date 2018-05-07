rm(list = ls()) #mainSsDbmh.R
library(RJafroc)

OptimisticScenario <-  FALSE
fileName <- "VanDyke.lrc" # fileName <- "Franken1.lrc"
cat("File name = ", fileName, "\n")
if (fileName == "VanDyke.lrc") dataset <- dataset02 else dataset <- dataset03

retDbm <- StSignificanceTesting(
  dataset02, 
  FOM = "Wilcoxon", 
  method = "DBMH") 
varYTR <- retDbm$varComp$varComp[3]
varYTC <- retDbm$varComp$varComp[4]
varYEps <- retDbm$varComp$varComp[6]
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
cat("observed effect size =", effectSize, "\n")
effectSize <- abs(effectSize)
sigma <- (retDbm$ciDiffTrtRRRC$`CI Upper`
          -retDbm$ciDiffTrtRRRC$`CI Lower`)/4
if (OptimisticScenario == TRUE) {
  if (fileName == "VanDyke.lrc") {
    effectSize <- effectSize -2*sigma
  }
  if (fileName == "Franken1.lrc") {
    effectSize <- effectSize +2*sigma
  }
}

cat("p-value = ", retDbm$pRRRC, 
    "\nanticipated effectSize = ", effectSize, 
    "\nCI Lower =", retDbm$ciDiffTrtRRRC$`CI Lower`,
    "\nCI Upper =", retDbm$ciDiffTrtRRRC$`CI Upper`, "\n")
powTab <- SsPowerTable(
  effectSize = effectSize, 
  method = "DBMH", 
  varYTR = varYTR, 
  varYTC = varYTC, 
  varYEps = varYEps)
print(powTab)