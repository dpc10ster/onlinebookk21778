rm( list = ls())#mainEmpVsFit.R
library(RJafroc);library(ggplot2)
seed <- 10;set.seed(seed)
mu <- 2; sigma <- 1.5
cat("Population AUC = ", 
    pnorm(mu/sqrt(1+sigma^2)), "\n")
K1 <- 500; K2 <- 500
fp <- rnorm(K1);tp <- rnorm(K2, mu, sigma)
zetas <- c(-Inf, 1.5, 2, 2.5, 3, 4, Inf)
fp1 <- as.numeric(cut(fp, zetas))
tp1 <- as.numeric(cut(tp, zetas))
rocData1 <- Df2RJafrocDataset(fp1, tp1)
plotEmp1 <- PlotEmpiricalOperatingCharacteristics(
  rocData1, 1, 1)
print(plotEmp1$Plot)
empAuc1 <- UtilFigureOfMerit(
  rocData1, FOM = "Wilcoxon")
cat("Emp. AUC bunched data = ", empAuc1, "\n")
Fit1 <- FitCbmRoc(rocData1)
print(Fit1$fittedPlot)
cat("CBM AUC, bunched data =" , Fit1$AUC,"\n")
zetas <- c(-Inf, -0.5, 0, 1, 1.5, 2, Inf)
fp2 <- as.numeric(cut(fp, zetas))
tp2 <- as.numeric(cut(tp, zetas))
rocData2 <- Df2RJafrocDataset(fp2, tp2)
plotEmp2 <- PlotEmpiricalOperatingCharacteristics(
  rocData2, 1, 1)
print(plotEmp2$Plot)
empAuc2 <- UtilFigureOfMerit(
  rocData2, FOM =
    "Wilcoxon")
cat("Emp. AUC well-spaced = ", empAuc2, "\n")
Fit2 <- FitCbmRoc(rocData2)
print(Fit2$fittedPlot)
cat("CBM AUC, well-spaced data =" , Fit2$AUC,"\n")