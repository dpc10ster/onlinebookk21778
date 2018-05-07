#mainPowerObservedEffectSize.R
# direct comparision of powers for observed effect sizes
rm(list = ls()) 
library(ggplot2)
library(RJafroc)

# included datasets
fileNames <-  c("TONY", "VD", "FR", 
                "FED", "JT", "MAG", 
                "OPT", "PEN", "NICO",
                "RUS", "DOB1", "DOB2", 
                "DOB3", "FZR")

f <- 4
fileName <- fileNames[f]
theData <- get(sprintf("dataset%02d", f))

retFileName <- paste0("allResults", fileName) 
sysAnalFileName <- system.file(
  "ANALYZED/RSM6", 
  retFileName, 
  package = "RJafroc")
load(sysAnalFileName)
I <- allResults[[1]]$I
J <- allResults[[1]]$J


# FED data has 5 modalities; we choose to analyze the first two
i1 <- 1;i2 <- 2 
cat("NH i1 = ", i1, "NH i2 = ", i2, "\n")
selectJ <- c(1, 2, 3, 4)
frocData <- DfExtractDataset(
  theData, trts = c(i1, i2), rdrs = selectJ)
J <- length(frocData$readerID)
K <- dim(frocData$NL)[3]

JTest <- 5;KTest <- 200
resROC <- StSignificanceTesting(
  frocData, 
  FOM = "HrAuc", 
  method = "DBMH", 
  option = "RRRC")
varCompROC <- resROC$varComp

effectSizeROC <- resROC$ciDiffTrtRRRC$Estimate

reswAFROC <- StSignificanceTesting(
  frocData, 
  FOM = "wAFROC", 
  method = "DBMH", 
  option = "RRRC")
varCompwAFROC <-reswAFROC$varComp 

effectSizewAFROC <- reswAFROC$ciDiffTrtRRRC$Estimate

cat("JTest = ", JTest, "KTest = ", KTest, "\n")
varYTR <- varCompROC$varComp[3]
varYTC <- varCompROC$varComp[4]
varYEps <- varCompROC$varComp[6]
powerROC <- SsPowerGivenJK(
  JTest, 
  KTest, 
  alpha = 0.05, 
  effectSize = effectSizeROC, 
  option = "RRRC",
  method = "DBMH", 
  varYTR = varYTR, 
  varYTC = varYTC, 
  varYEps = varYEps)$powerRRRC

varYTR <- varCompwAFROC$varComp[3]
varYTC <- varCompwAFROC$varComp[4]
varYEps <- varCompwAFROC$varComp[6]
powerwAFROC <- SsPowerGivenJK(
  JTest, KTest, 
  alpha = 0.05, 
  effectSize = effectSizewAFROC, 
  option = "RRRC",
  method = "DBMH", 
  varYTR = varYTR, 
  varYTC = varYTC, 
  varYEps = varYEps)$powerRRRC

cat("ROC effect-size = ,", effectSizeROC, 
    "wAFROC effect-size = ,", effectSizewAFROC, 
    "Power ROC, wAFROC:", powerROC, ",", powerwAFROC, "\n")

