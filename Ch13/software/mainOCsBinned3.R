# mainOCsBinned.R
rm(list = ls())
library(RJafroc);library(ggplot2)

seed <- 1;set.seed(seed)
mu <- 1;nu <- 1;lambda <- 2;zeta1 <- -1;K1 <- 5;K2 <- 7 # Fig. 13.5 (a-c)
mu <- 1;nu <- 1;lambda <- 2;zeta1 <- -1;K1 <- 50;K2 <- 70 # Fig. 13.5 (d-f)
# mu <- 0.1;lambda <- 1;nu <- 1 ;zeta1 <- -1;K1 <- 50;K2 <- 70 # Fig. 13.1
# mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- -1;K1 <- 5;K2 <- 7 # Column 1 in Fig. 13.2
# mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- -1;K1 <- 50;K2 <- 70 # Column 2 in Fig. 13.2
# mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- +1;K1 <- 50;K2 <- 70 # Column 3 in Fig. 13.2
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1))

frocDataRaw <- SimulateFrocDataset(
  mu = mu, lambda = lambda, nu = nu, 
  I = 1, J = 1, K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1
)
frocDataBin <- DfBinDataset(frocDataRaw, opChType = "FROC")

plotFROC <- PlotEmpiricalOperatingCharacteristics(frocDataBin, trts= 1, rdrs = 1, opChType = "FROC")
plotFROC1 <- plotFROC$Plot + scale_color_manual(values = "black")  +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
plotFROC1$layers[[1]]$aes_params$size <- 2 # line
plotFROC1$layers[[2]]$aes_params$size <- 5 # points
print(plotFROC1)

afrocDataBin <- DfBinDataset(frocDataRaw, opChType = "AFROC")
plotAFROC <- PlotEmpiricalOperatingCharacteristics(afrocDataBin, trts= 1, rdrs = 1, opChType = "AFROC")
plotAFROC1 <- plotAFROC$Plot + scale_color_manual(values = "black")  +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
plotAFROC1$layers[[1]]$aes_params$size <- 2 # line
plotAFROC1$layers[[2]]$aes_params$size <- 5 # points
print(plotAFROC1)

rocDataBin <- DfBinDataset(frocDataRaw, opChType = "ROC")
plotROC <- PlotEmpiricalOperatingCharacteristics(rocDataBin, trts= 1, rdrs = 1, opChType = "ROC")
plotROC1 <- plotROC$Plot + scale_color_manual(values = "black")  +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
plotROC1$layers[[1]]$aes_params$size <- 2 # line
plotROC1$layers[[2]]$aes_params$size <- 5 # points
print(plotROC1)
