# mainQuantifySearchPerformance.R
rm(list = ls())
library(RJafroc)
library(ggplot2)
library(grid)
source("AucsRSM.R")

a <- 2.5; b <- 1
mu <- 2; lambda <- 1; nu <- 1; 
nLesDistr <- c(1, 1); if (length(nLesDistr) == 2) dim(nLesDistr) <- c(1, 2) 

plotZeta <- seq(-20, 20, by = 0.01)
lambdaP <- lambda / mu
if (abs(nu * mu) <= 1e-6 ) nuP <- 1e-6 else nuP <- (1-exp(-nu * mu))

FPFBM <- 1 - pnorm(plotZeta)
TPFBM <- 1 - pnorm(plotZeta, mean = a/b, sd = 1/b)
plotBM <- data.frame(FPF = FPFBM, TPF = TPFBM)

FPFRSM <- sapply(plotZeta, xROC, lambdaP = lambdaP)
TPFRSM <- sapply(plotZeta, yROC, mu = mu, lambdaP = lambdaP, 
                nuP = nuP, lesionDistribution = nLesDistr)
plotRSM <- data.frame(FPF = FPFRSM, TPF = TPFRSM)
dashedRSM <- data.frame(FPF = c(FPFRSM[1], 1), TPF = c(TPFRSM[1], 1))

fpfMax <- max(FPFRSM)
tpfMax <- max(TPFRSM)
if (fpfMax < 0.99){
  fpfCross <- (fpfMax + tpfMax) / 2
  tpfCross <- fpfCross
  endPointRSM <- data.frame(FPF = fpfMax, TPF = tpfMax)
  endPointBM <- data.frame(FPF = 1, TPF = 1)
  ds <- data.frame(FPF = c(fpfMax, fpfCross), TPF = c(tpfMax, tpfCross))
  diagonal <- data.frame(FPF = c(0, 1), TPF = c(0, 1))
  aText <- data.frame(FPF = 0.1, TPF = 0.95)
  bText <- data.frame(FPF = 0.2, TPF = 0.8)
  cText <- data.frame(FPF = 0.5, TPF = 0.45)
  dsText <- data.frame(FPF = (fpfMax + fpfCross)/2 + 0.05, TPF = (tpfMax + tpfCross)/2)
  compPlot <- ggplot() + 
    geom_line(mapping = aes(x = FPF, y = TPF), data = plotRSM, size = 2) + 
    geom_line(mapping = aes(x = FPF, y = TPF), data = plotBM, size = 2) + 
    geom_line(data = dashedRSM, aes(x = FPF, y = TPF), linetype = 3, size = 2) +
    geom_point(data = endPointRSM, mapping = aes(x = FPF, y = TPF), shape = 15, size = 7) +
    geom_point(data = endPointBM, mapping = aes(x = FPF, y = TPF), shape = 16, size = 7) +
    geom_line(data = ds, mapping = aes(x = FPF, y = TPF), linetype = 2) + 
    geom_line(data = diagonal, mapping = aes(x = FPF, y = TPF), linetype = 2) + 
    geom_text(data = dsText, mapping = aes(x = FPF, y = TPF), label = "d[s]", parse = TRUE, size = 20) +
    geom_text(data = aText, mapping = aes(x = FPF, y = TPF), label = "a", parse = TRUE, size = 20) +
    geom_text(data = bText, mapping = aes(x = FPF, y = TPF), label = "b", parse = TRUE, size = 20) +
    geom_text(data = cText, mapping = aes(x = FPF, y = TPF), label = "c", parse = TRUE, size = 20) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    coord_fixed(ratio = 1) #
  compPlot <- compPlot + theme(axis.title.y = element_text(size = 25,face="bold"),
                             axis.title.x = element_text(size = 30,face="bold"))
  compPlot <- ggplotGrob(compPlot)
  compPlot$layout$clip[compPlot$layout$name=="panel"] <- "off"
  grid.draw(compPlot)
}


