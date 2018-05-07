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
  afrocPlot <- ret$AFROCPlot + 
    scale_color_manual(values = "black") + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold")) 
  afrocPlot$layers[[1]]$aes_params$size <- 2 # line
  afrocPlot$layers[[2]]$aes_params$size <- 2 # line
  print(afrocPlot) # AFROC plot
  cat("mu = ", mu,
      "\nlambda = ", lambda,
      "\nnu = ", nu, 
      "\nAUC = ", ret$aucAFROC,
      "\nfpfMax = ", max(ret$AFROCPlot$data$FPF),
      "\nllfMax = ", max(ret$AFROCPlot$data$LLF),"\n")
}