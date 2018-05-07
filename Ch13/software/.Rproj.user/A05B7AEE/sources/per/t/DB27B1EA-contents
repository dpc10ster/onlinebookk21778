rm(list = ls()) # mainOCsRaw.R
library(RJafroc);library(ggplot2)

seed <- 1;set.seed(seed)
mu <- 1;lambda <- 1;nu <- 1
zeta1 <- -1;K1 <- 5;K2 <- 7
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1))

frocDataRaw <- SimulateFrocDataset(
  mu = mu, lambda = lambda, nu = nu, I = 1, J = 1,
  K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1)

plotFROC <- PlotEmpiricalOperatingCharacteristics(
  dataset = frocDataRaw,
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

plotAFROC <- PlotEmpiricalOperatingCharacteristics(
  dataset = frocDataRaw, 
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

plotROC <- PlotEmpiricalOperatingCharacteristics(
  dataset = frocDataRaw,
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

retRocRaw <- DfFroc2Roc(frocDataRaw)
