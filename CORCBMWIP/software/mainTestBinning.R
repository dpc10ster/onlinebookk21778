rm(list = ls())
library(RJafroc)
library(ggplot2)
library(bbmle)
options(digits = 6)

#dataset <- dataset05;dataset <- DfFroc2Roc(dataset)
seed <- 123
dataset <- SimulateRocDataset(K1 = 5000, K2 = 7000, a = 1, b = 0.5, seed = 123)
fomOrg <- as.matrix(UtilFigureOfMerit(dataset, fom = "Wilcoxon"), nrow = 2, ncol = 9)
print(fomOrg)
cat("mean, sd = ", mean(fomOrg), sd(fomOrg), "\n")
I <- length(dataset$NL[,1,1,1])
J <- length(dataset$NL[1,,1,1])
K1 <- length(dataset$NL[1,1,,1])
K2 <- length(dataset$LL[1,1,,1])
K1 <- K1 - K2
desiredNumBins <- 7
numZeta <- desiredNumBins - 1
localMax <- array(-1,dim = c(I,J))
sSave <- array(dim = c(I,J))
zetasArr <- array(dim = c(I,J,numZeta))
datasetB <- dataset
for (i in 1:I)
{
  for (j in 1:J)
  { 
    numZeta0 <- numZeta
    NL <- dataset$NL[i,j,1:K1,1]
    LL <- dataset$LL[i,j,1:K2,1]
    nLlL <- c(NL,LL)
    # need to remove lowest value, as this gives (1,1) point
    candidateZetas <-  sort(unique(nLlL))[-1]
    el <- length(candidateZetas)
    if (el < numZeta) {
      numZeta0 <- el - 1
      sample <- combn(candidateZetas, numZeta0)
    } else {
      # if more than 20 candidates, need to trim
      if (el > 20) {
        byDivisor <- 10
        while (1) {
          by <- as.integer(el/byDivisor)
          candidateZetasTrim <- candidateZetas[seq(from = 1, to = el, by = by)]
          sample <- combn(candidateZetasTrim, numZeta)
          if (length(sample[1,]) > 200) {
            byDivisor <- byDivisor - 1
          } else break
        }
      } else sample <- combn(candidateZetas, numZeta)
    }
    for (s in 1:length(sample[1,]))
    {
      zetas1 <- sort(sample[,s])
      zetas <- c(-Inf,zetas1,+Inf)
      nLlLB <- cut(nLlL, zetas, labels = FALSE, right = FALSE)
      datasetB$NL[i,j,1:K1,1] <- nLlLB[1:K1]
      datasetB$LL[i,j,1:K2,1] <- nLlLB[(K1+1):(K1+K2)]
      fom <- UtilFigureOfMerit(datasetB, fom = "Wilcoxon")[i,j]
      if (fom > localMax[i,j]){
        sSave[i,j] <- s
        localMax[i,j] <- fom
        zetasArr[i,j,] <- zetas1
      }
    }
    next
  }
}

datasetB <- dataset
for (i in 1:I)
{
  for (j in 1:J)
  { 
    NL <- dataset$NL[i,j,1:K1,1]
    LL <- dataset$LL[i,j,1:K2,1]
    zetas <- c(-Inf,zetasArr[i,j,],Inf)
    nLlLB <- cut(nLlL, zetas, labels = FALSE, right = FALSE)
    datasetB$NL[i,j,1:K1,1] <- nLlLB[1:K1]
    datasetB$LL[i,j,1:K2,1] <- nLlLB[(K1+1):(K1+K2)]
  }
}
fom <- UtilFigureOfMerit(datasetB, fom = "Wilcoxon")
print(fom)
cat("mean, sd = ", mean(fom), sd(fom), "\n")
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