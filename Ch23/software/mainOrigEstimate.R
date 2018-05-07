rm(list = ls())
library(RJafroc)
library(doParallel)
library(doRNG)
library(foreach)
library(mvtnorm)
library(lmf)

stop("fix or delete me")
fileName <- "NICO"
rocData <- loadDataFile(fileName, "./")
retOrh <- UtilVarianceComponents(rocData, "Wilcoxon", method = "ORH")
J <- length(rocData$readerID)
K <- dim(rocData$NL)[3]
K2 <- dim(rocData$LL)[3]
K1 <- K - K2

seed <- 1
cl <- makeCluster(detectCores())
registerDoParallel(cl)
B <- 2000
retBs <- foreach(b  = 1:B, .combine = "rbind", .packages = "RJafroc", .options.RNG = seed) %dorng% {
  jBs <- ceiling(runif(J, 0, J))
  k1Bs <- ceiling(runif(K1, 0, K1))
  k2Bs <- ceiling(runif(K2, 0, K2))
  rocDataTemp <- rocData
  rocDataTemp$NL <- rocData$NL[ , jBs, c(k1Bs, k2Bs + K1), , drop = FALSE]
  rocDataTemp$LL <- rocData$LL[ , jBs, k2Bs, , drop = FALSE]
  ret <- UtilVarianceComponents(rocDataTemp, "Wilcoxon", method = "ORH")
  # c(ret$cov2, ret$varEps, jBs, k1Bs, k2Bs)
  c(ret$cov2, ret$varEps)
}
stopCluster(cl)
cat("cov2 and var of Orig. data are:", retOrh$cov2, retOrh$varEps, "\n")
cat("95% CI of cov2 is: (", quantile(retBs[ , 1], 0.025), ",", quantile(retBs[ , 1], 0.975), ")\n")
cat("95% CI of var is: (", quantile(retBs[ , 2], 0.025), ",", quantile(retBs[ , 2], 0.975), ")\n")
