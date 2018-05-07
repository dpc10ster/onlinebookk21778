rm(list = ls()) # mainRsmPlotsAfroc.R
library(RJafroc)
muArr <- c(0.001,1,2,3,4,5)
lambda <- 1;nu <- 1; Lmax <- 10 
for (i in 1:length(muArr)) {
  mu <- muArr[i]
  ret <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, 
    type = "AFROC", 
    lesionDistribution = Lmax) 
  afrocPlot <- ret$AFROCPlot
  print(afrocPlot)
  cat("mu = ", mu,
      "\nlambda = ", lambda,
      "\nnu = ", nu, 
      "\nAUC = ", ret$aucAFROC,
      "\nfpfMax = ", max(ret$AFROCPlot$data$FPF),
      "\nllfMax = ", max(ret$AFROCPlot$data$LLF),"\n")
}