rm(list = ls()) # mainRsmPlotsFROC.R
library(RJafroc)
muArr <- c(0.001,1,2,3,4,5)
lambda <- 1;nu <- 1; Lmax <- 1 
for (i in 1:length(muArr)) {
  mu <- muArr[i]
  ret <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, 
    type = "FROC", lesionDistribution = Lmax,
    llfRange = c(0,1)) 
  frocPlot <- ret$FROCPlot + 
    scale_color_manual(values = "black") + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold")) 
  frocPlot$layers[[1]]$aes_params$size <- 2 # line
  print(frocPlot) 
  cat("mu = ", mu,
      "\nlambda = ", lambda,
      "\nnu = ", nu, 
      "\nAUC = ", ret$aucFROC,
      "\nnlfMax = ", max(ret$FROCPlot$data$NLF),
      "\nllfMax = ", max(ret$FROCPlot$data$LLF),"\n")
}
