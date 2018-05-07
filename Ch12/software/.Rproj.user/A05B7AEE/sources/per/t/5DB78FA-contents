# MainFrocCurvePop.R
rm(list = ls())
library(RJafroc);library(ggplot2)

seed <- 1;set.seed(seed)
# 0.5, 1, 2, 10, 0.1
mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- -Inf
K1 <- 10000;K2 <- 10000
K1 <- 4;K2 <- 3
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1))

frocDataRaw <- SimulateFrocDataset(
  mu = mu, lambda = lambda, nu = nu, 
  I = 1, J = 1, K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1
)

plotFROC <- PlotEmpiricalOperatingCharacteristics(
  dataset = frocDataRaw, trts= 1, rdrs = 1, opChType = "FROC")
p <- plotFROC$Plot +
  scale_color_manual(values = "black")  +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
p$layers[[1]]$aes_params$size <- 2 # line
p$layers[[2]]$aes_params$size <- 5 # points
print(p)
