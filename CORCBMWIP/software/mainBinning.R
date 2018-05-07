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

AddArguments <- function(f,n){
  # adds n arguments to a function f; returns that new function 
  t <- paste("arg <- alist(", paste(names(formals(f)), collapse="=, "), "=, ",
             paste(sapply(1:n, function(i) paste("zetaFwd", i, "=", sep = "")), collapse = ", "),
             ")", sep= "")
  formals(f) <- eval(parse(text = t))
  return(f)
}

#data = list(NL = NL, LL = LL, fomOrg = fomOrg, K1 = K1, K2 = K2)
minimizationFunction <- function (NL, LL, fomOrg, K1, K2){
  allParameters <- names(formals())
  zetaPos <- regexpr("zeta", allParameters)
  zetaFwd <- unlist(mget(allParameters[which(zetaPos == 1)]))
  zetas <- InverseZetas(zetaFwd)
  #cat("zetas = ", zetas, "\n")
  zetas <- c(-Inf, zetas, Inf)

  nl <- cut(NL, zetas, labels = FALSE, right = FALSE)
  ll <- cut(LL, zetas, labels = FALSE, right = FALSE)
  dataset <- Df2RJafrocDataset(nl, ll)
  fomBin <- UtilFigureOfMerit(dataset, fom = "Wilcoxon")
  diff <- abs(fomOrg - fomBin)
  diff <- fomOrg - fomBin
  #cat("diff = ", diff,"\n")
  return(diff)
}


seed <- 123
a <- 1.5;b <- 0.5
dataset <- SimulateRocDataset(K1 = 5000, K2 = 7000, a = a, b = b, seed = 123)
dataset <- dataset05
K1 <- length(dataset$NL[1,1,,1])
K2 <- length(dataset$LL[1,1,,1])
K1 <- K1 - K2
NL <- dataset$NL[1,1,1:K1,1]
LL <- dataset$LL[1,1,1:K2,1]
fomOrg <- UtilFigureOfMerit(dataset, fom = "Wilcoxon")
cat("fomOrg = ", fomOrg,"\n")
desiredNumBins <- 5
zetasIni <- seq(-2, a/b+2, length.out = 7) + 1
cat("zetasIni = ", zetasIni, "\n")
zetasFwd <- ForwardZetas(zetasIni)
namesVector <- c(paste0("zetaFwd", 1:length(zetasIni)))
parameters <- as.list(zetasFwd)
names(parameters) <- namesVector
minimizationFunctionNew <- AddArguments(minimizationFunction, length(zetasFwd))

ret <- mle2(minimizationFunctionNew, start = parameters, 
            data = list(NL = NL, LL = LL, fomOrg = fomOrg, K1 = K1, K2 = K2))
cat("value =", InverseZetas(ret@coef),"\n")
cat("value =", ret@min,"\n")