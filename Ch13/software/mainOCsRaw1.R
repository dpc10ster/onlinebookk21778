# mainOCsRaw.R
rm(list = ls())
library(RJafroc);library(ggplot2)

seed <- 1;set.seed(seed)
# mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- +1;K1 <- 50;K2 <- 70
# mu <- 0.1;lambda <- 1;nu <- 1 ;zeta1 <- -1;K1 <- 50;K2 <- 70 # Table 13.1
mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- -1;K1 <- 5;K2 <- 7 # Column 1 in Table 13.2
mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- -1;K1 <- 50;K2 <- 70 # Column 2 in Table 13.2
#mu <- 1;lambda <- 1;nu <- 1 ;zeta1 <- +1;K1 <- 50;K2 <- 70 # Column 3 in Table 13.2
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1))

frocDataRaw <- SimulateFrocDataset(
  mu = mu, lambda = lambda, nu = nu, 
  I = 1, J = 1, K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1
)

# p <- PlotEmpiricalOperatingCharacteristics(dataset = frocDataRaw, trts= 1, rdrs = 1, opChType = "FROC")
# p <- p$Plot + scale_color_manual(values = "black")  +
#   theme(legend.position="none") +
#   theme(axis.title.y = element_text(size = 25,face="bold"),
#         axis.title.x = element_text(size = 30,face="bold"))
# #plotFROC1$layers[[1]]$aes_params$size <- 4 # line
# p$layers[[1]]$aes_params$size <- 2 # line
# p$layers[[2]]$aes_params$size <- 5 # points
# print(p)

p <- PlotEmpiricalOperatingCharacteristics(dataset = frocDataRaw, trts= 1, rdrs = 1, opChType = "AFROC")
p <- p$Plot + scale_color_manual(values = "black")  +
  theme(legend.position="none") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
p$layers[[1]]$aes_params$size <- 2 # line
p$layers[[2]]$aes_params$size <- 5 # points
print(p) 

# p <- PlotEmpiricalOperatingCharacteristics(dataset = frocDataRaw, trts= 1, rdrs = 1, opChType = "ROC")
# p <- p$ROCPlot + scale_color_manual(values = "black")  +
#   theme(legend.position="none") +
#   theme(axis.title.y = element_text(size = 25,face="bold"),
#         axis.title.x = element_text(size = 30,face="bold"))
# p$layers[[1]]$aes_params$size <- 2 # line
# p$layers[[2]]$aes_params$size <- 5 # points
# print(p) 
# 
# retRocRaw <- DfFroc2HrRoc(frocDataRaw)
