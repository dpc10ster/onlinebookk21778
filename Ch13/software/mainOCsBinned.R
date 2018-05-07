rm(list = ls()) # mainOCsBinned.R
library(RJafroc);library(ggplot2)

seed <- 1;set.seed(seed)
mu <- 1;lambda <- 1;nu <- 1
zeta1 <- -1;K1 <- 50;K2 <- 70
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1))

frocDataRaw <- SimulateFrocDataset(
  mu = mu, lambda = lambda, nu = nu, I = 1, J = 1,
  K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1)

frocDataBin <- DfBinDataset(frocDataRaw, desiredNumBins = 5, opChType = "FROC")
plotFROC <- PlotEmpiricalOperatingCharacteristics(
  dataset = frocDataBin,
  trts= 1,
  rdrs = 1,
  opChType = "FROC")
p <- plotFROC$Plot +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
p$layers[[1]]$aes_params$size <- 2 # line
p$layers[[2]]$aes_params$size <- 5 # points
print(p)

afrocDataRaw <- DfFroc2Afroc(frocDataRaw)
afrocDataBin <- DfBinDataset(afrocDataRaw, desiredNumBins = 5, opChType = "AFROC")
plotAFROC <- PlotEmpiricalOperatingCharacteristics(
  dataset = afrocDataBin, 
  trts= 1, 
  rdrs = 1, 
  opChType = "AFROC" 
)
p <- plotAFROC$Plot +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
p$layers[[1]]$aes_params$size <- 2 # line
p$layers[[2]]$aes_params$size <- 5 # points
print(p)

rocDataRaw <- DfFroc2Roc(frocDataRaw)
rocDataBin <- DfBinDataset(rocDataRaw, desiredNumBins = 5, opChType = "ROC")
plotROC <- PlotEmpiricalOperatingCharacteristics(
  dataset = rocDataBin,
  trts= 1,
  rdrs = 1,
  opChType = "ROC")
p <- plotROC$Plot +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
p$layers[[1]]$aes_params$size <- 2 # line
p$layers[[2]]$aes_params$size <- 5 # points
print(p)
