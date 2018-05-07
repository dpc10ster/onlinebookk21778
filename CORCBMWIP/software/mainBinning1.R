rm(list = ls())
library(RJafroc)
library(ggplot2)
library(bbmle)
options(digits = 6)
ForwardZetas <- function(zetas) {
  .Call('RJafroc_ForwardZetas', PACKAGE = 'RJafroc', zetas)
}

InverseZetas <- function(zetasFwd) {
  .Call('RJafroc_InverseZetas', PACKAGE = 'RJafroc', zetasFwd)
}

minimizationFunction <- function (zetaFwd, NL, LL, fomOrg, K1, K2){
  zetas <- InverseZetas(zetaFwd)
  #cat("zetas = ", zetas, "\n")
  zetas <- c(-Inf, zetas, Inf)
  #if (max(zetas) > max(c(NL, LL))) return (1000)
  nl <- cut(NL, zetas, labels = FALSE, right = FALSE)
  ll <- cut(LL, zetas, labels = FALSE, right = FALSE)
  dataset <- Df2RJafrocDataset(nl, ll)
  fomBin <- UtilFigureOfMerit(dataset, fom = "Wilcoxon")
  #diff <- abs(fomOrg - fomBin)
  #diff <- fomOrg - fomBin
  #cat("diff = ", diff,"\n")
  return(-fomBin)
}


#seed <- 123
#dataset <- SimulateRocDataset(K1 = 50, K2 = 70, a = 1, b = 0.5, seed = 123)
dataset <- dataset05;dataset <- DfFroc2Roc(dataset)
fomOrg <- UtilFigureOfMerit(dataset, fom = "Wilcoxon")
K1 <- length(dataset$NL[1,1,,1])
K2 <- length(dataset$LL[1,1,,1])
K1 <- K1 - K2
datasetB <- dataset
for (i in 1:2)
{
  for (j in 1:9)
  {    
    cat("fomOrg[i,j] = ", fomOrg[i,j],"\n")
    NL <- dataset$NL[i,j,1:K1,1]
    LL <- dataset$LL[i,j,1:K2,1]
    cat("unique = ", sort(unique(c(NL, LL))), "\n")
    nLlL <- c(NL,LL)
    #zetasIni <- c(min(nLlL) - 1,  seq(min(nLlL), max(nLlL), length.out = 6), max(nLlL)+ 1)
    zetasIni <- seq(min(nLlL) + 0.1, max(nLlL) - 0.1, length.out = 6)
    cat("zetasIni = ", zetasIni, "\n")
    zetasFwd <- ForwardZetas(zetasIni)
    ret <- optim(zetasFwd, minimizationFunction, gr = NULL, NL, LL, fomOrg, K1, K2, method = "BFGS",
                 control = list(parscale=rep(0.9,length(zetasIni))))
    cat(InverseZetas(ret$par), "\n", -ret$value, "\n")

    nLlLB <- cut(nLlL, zetas, labels = FALSE, right = FALSE)
    datasetB$NL[i,j,1:K1,1] <- nLlLB[1:K1]
    datasetB$NL[i,j,1:K2,1] <- nLlLB[(K1+1):(K1+K2)]
    next
  }
}
fomBin <- UtilFigureOfMerit(dataset, fom = "Wilcoxon")
cat("fomBin = ", fomBin,"\n")
cat("fomOrg - fomBin = ", fomOrg - fomBin,"\n")
zetasIni <- seq(-2, 6, length.out = 6) 
cat("zetasIni = ", zetasIni, "\n")
zetasFwd <- ForwardZetas(zetasIni)
ret <- optim(zetasFwd, minimizationFunction, gr = NULL, NL, LL, fomOrg, K1, K2, method = "SANN",
             control = list(parscale=rep(0.01,length(zetasIni))))
cat(InverseZetas(ret$par), "\n", ret$value, "\n")

