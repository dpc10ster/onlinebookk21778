rm(list = ls())
library(RJafroc)
library(ggplot2)
#source("DfBinDataset3.R")
options(digits = 6)

# JT data
dataset <- dataset05;dataset <- DfFroc2Roc(dataset)
# seed <- 123
# dataset <- SimulateRocDataset(K1 = 5000, K2 = 7000, a = 1, b = 0.5, seed = 123)
I <- length(dataset$NL[,1,1,1])
J <- length(dataset$NL[1,,1,1])

datasetB <- DfBinDataset(dataset)

if ((I == 1) && (J == 1)) {
  x <- PlotEmpiricalOperatingCharacteristics(dataset)$Plot
  y <- PlotEmpiricalOperatingCharacteristics(datasetB)$Points
  fpf <- y$genAbscissa[-1];fpf <- fpf[-length(fpf)]
  tpf <- y$genOrdinate[-1];tpf <- tpf[-length(tpf)]
  plotOpPnts <- rbind(data.frame(fpf = fpf, tpf = tpf))
  x <- x + geom_point(data = plotOpPnts, aes(x = fpf, y = tpf), size = 4)
  print(x)
} else {
  x <- PlotEmpiricalOperatingCharacteristics(dataset, trts = c(1,2), rdrs = seq(1,9))$Plot;print(x)
  x <- PlotEmpiricalOperatingCharacteristics(datasetB, trts = c(1,2), rdrs = seq(1,9))$Plot;print(x)
}