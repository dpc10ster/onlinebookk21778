rm(list = ls()) # mainNoMarks.R
library(RJafroc)

Lmax <- 1;mu <- 0.001
lambda <- 0.000001;nu <- 1

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
#print(retFroc$FROCPlot) # FROC plot
afrocPlot <- retAfroc$AFROCPlot + scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) 
afrocPlot$layers[[1]]$aes_params$size <- 2 # line
afrocPlot$layers[[2]]$aes_params$size <- 2 # line
print(afrocPlot) # AFROC plot
print(afrocPlot) # AFROC plot
#print(retRoc$ROCPlot) # ROC plot
cat("mu = ", mu,
    "\nlambda = ", lambda,
    "\nnu = ", nu, 
    "\naucFroc = ", retFroc$aucFROC,
    "\naucRoc = ", retRoc$aucROC,
    "\naucAfroc = ", retAfroc$aucAFROC,"\n")
