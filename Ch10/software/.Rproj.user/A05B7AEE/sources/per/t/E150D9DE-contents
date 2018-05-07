rm(list = ls()) # mainSingleTreatment.R
library(RJafroc)
source("Wilcoxon.R")
source("CovJk.R")
source("CovDL.R")
source("SingleTreatmentAnalysis.R")
alpha <- 0.05
ROC <- FALSE
if (ROC) {
  #fileName <- "Franken1.lrc"
  fileName <- "VanDyke.lrc"
  rocData <- DfReadDataFile(
    fileName, 
    format = "MRMC")
} else {
  fileName <- "CXRinvisible3-20mm.xlsx"
  frocData <- DfReadDataFile(
    fileName, 
    format = "JAFROC")
  rocData <- DfFroc2Roc(frocData)
}
cat("data file = ", fileName, "\n")

i <- 1 # select the treatment to be analyzed
# extract the first treatment
rocData1T <- 
  DfExtractDataset(rocData, trts = i)
fomArray <- 
  UtilFigureOfMerit(
    rocData1T, FOM = "Wilcoxon")
thetaDot <- mean(fomArray[i, ])
mu0 <- 0.583422;mu0 <- 0.6#mu0 <- 0.583422#
ret <- SingleTreatmentAnalysis(
  rocData1T, 
  mu0, 
  covEstMthd = "Jackknife", 
  alpha = alpha)
cat("The NH is that thetaDot = mu0, where thetaDot= ", 
    thetaDot, "\nand mu0 = ", mu0,"\n")
cat("The mean FOM for the anal2zed treatment is:", thetaDot,"\n")
cat("The", 100 * (1 - alpha), "% CI for the preceding value is:", "(", ret$ci[1], ",", ret$ci[2], ")\n")
cat("The t-statistic to test\nH0: (analyzed treatment = standard) is:", ret$tStat, 
    "\nand the and p-value is ", ret$pVal, "\n")
cat("The difference in rdr.avg  minus standard = ",thetaDot - mu0, "\n")
cat("The", 100 * (1 - alpha), "% CI of the preceding value is", "(", ret$ciDiff[1], ",", ret$ciDiff[2], ")\n")
