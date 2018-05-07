# mainNoMarks.R
rm(list = ls())
library(RJafroc)

Lmax <- 1;mu <- 0.001;lambda <- 0.000001;nu <- 1

retFroc <- PlotRsmOperatingCharacteristics(
  mu, lambda, nu, 
  type = "FROC", lesionDistribution = Lmax, 
  llfRange = c(0,1)
) 
retRoc <- PlotRsmOperatingCharacteristics(
  mu, lambda, nu, 
  type = "ROC", lesionDistribution = Lmax
) 
retAfroc <- PlotRsmOperatingCharacteristics(
  mu, lambda, nu, 
  type = "AFROC", lesionDistribution = Lmax
) 

afrocPlot <- retAfroc$AFROCPlot
print(afrocPlot) # AFROC plot
cat("mu = ", mu,
    "\nlambda = ", lambda,
    "\nnu = ", nu, 
    "\naucFroc = ", retFroc$aucFROC,
    "\naucRoc = ", retRoc$aucROC,
    "\naucAfroc = ", retAfroc$aucAFROC,"\n")
