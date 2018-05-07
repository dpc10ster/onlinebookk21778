# mainFrocCurveBinned.R
rm(list = ls())
library(RJafroc);library(ggplot2)

muArr <-  c(0.5, 1, 2, 10, 0.01)
nbins <- c(6, 6, 6)
for (i in 1:3){
  seed <- 1;set.seed(seed)
  mu <- muArr[i];lambda <- 1;nu <- 1 ;zeta1 <- -1;K1 <- 50;K2 <- 70
  Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1))
  frocDataRaw <- SimulateFrocDataset(
    mu = mu, lambda = lambda, nu = nu, 
    I = 1, J = 1, K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1)
  
  frocDataBinned <- DfBinDataset(
    frocDataRaw, 
    desiredNumBins = nbins[i], 
    opChType = "FROC")
  
  plotFROC <- PlotEmpiricalOperatingCharacteristics(
    frocDataBinned, 
    trts= 1, 
    rdrs = 1, 
    opChType = "FROC")
  
  p <- plotFROC$Plot +   scale_color_manual(values = "black")  +
    theme(legend.position="none") +
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold")) +
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"))
  p$layers[[2]]$aes_params$size <- 5 # points
  print(p)
}