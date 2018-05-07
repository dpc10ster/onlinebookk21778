rm(list = ls())
library(RJafroc)
library(ggplot2)
library(bbmle)
options(digits = 6)

dataset <- SimulateRocDataset(K1 = 5000, K2 = 7000, a = 1, b = 0.5, seed = 123)
datasetB <- DfBinDataset(dataset, desiredNumBins = 7)
fomOrg <- as.matrix(UtilFigureOfMerit(dataset, fom = "Wilcoxon"), nrow = 2, ncol = 9)
print(fomOrg)
fomBinned <- as.matrix(UtilFigureOfMerit(datasetB, fom = "Wilcoxon"), nrow = 2, ncol = 9)
print(fomOrg)
cat("fomOrg = ", mean(fomOrg), "\n")
cat("fomBinned = ", mean(fomBinned), "\n")
x <- PlotEmpiricalOperatingCharacteristics(dataset)$Plot
y <- PlotEmpiricalOperatingCharacteristics(datasetB)$Points
fpf <- y$genAbscissa[-1];fpf <- fpf[-length(fpf)]
tpf <- y$genOrdinate[-1];tpf <- tpf[-length(tpf)]
plotOpPnts <- rbind(data.frame(fpf = fpf, tpf = tpf))
x <- x + geom_point(data = plotOpPnts, aes(x = fpf, y = tpf), size = 4)
print(x)
xx <- PlotEmpiricalOperatingCharacteristics(datasetB)
print(xx$Points)