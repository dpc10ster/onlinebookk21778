rm(list = ls())
library(mvtnorm)
library(RJafroc)
library(ggplot2)
library(binom)
library(doParallel)
library(doRNG)
library(foreach)

S <- 20
muX <- 1.5;muY <- 3
alphaX <- 0.4;alphaY <- 0.7
rhoNor <- 0.3;rhoAbn1 <- rhoNor;rhoAbn2 <- 0.8
rhoAbn3 <- mean(c(rhoAbn1, rhoAbn2))
sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
sigmaAbn1 <- rbind(c(1, rhoAbn1), c(rhoAbn1, 1))
sigmaAbn2 <- rbind(c(1, rhoAbn2), c(rhoAbn2, 1))
sigmaAbn3 <- rbind(c(1, rhoAbn3), c(rhoAbn3, 1))
# 50/50; 100/100; 1000/1000; 5000/5000
K1 <- 50;K2 <- 50
p200 <- (1 - alphaX) * (1 - alphaY)
p2X0 <- alphaX * (1 - alphaY)
p20Y <- (1 - alphaX) * alphaY
p2XY <- alphaX * alphaY

seed <- 123
cl <- makeCluster(detectCores())
registerDoParallel(cl)
ret <- foreach(s = 1:S, .combine = "rbind", .options.RNG = seed, .packages = c("mvtnorm", "binom", "RJafroc")) %dorng%{
  K2Sample <- sample(c("00", "X0", "0Y", "XY"), size = K2, replace = TRUE, prob = c(p200, p2X0, p20Y, p2XY))
  K200 <- sum(K2Sample == "00")
  K2X0 <- sum(K2Sample == "X0")
  K20Y <- sum(K2Sample == "0Y")
  K2XY <- sum(K2Sample == "XY")
  
  zk1 <- t(rmvnorm(K1, sigma = sigmaNor))
  zk200 <- t(rmvnorm(K200, mean = c(0, 0), sigma = sigmaAbn1))
  zk2X0 <- t(rmvnorm(K2X0, mean = c(muX, 0), sigma = sigmaAbn3))
  zk20Y <- t(rmvnorm(K20Y, mean = c(0, muY), sigma = sigmaAbn3))
  zk2XY <- t(rmvnorm(K2XY, mean = c(muX, muY), sigma = sigmaAbn2))
  
  zk2 <- cbind(zk200, zk2X0, zk20Y, zk2XY)
  
  # gridData <- expand.grid(x = seq(min(zk2[1, ]) - 0.5, max(zk2[1, ]) + 0.5, length.out = 200), y = seq(min(zk2[2, ]) - 0.5, max(zk2[2, ]), length.out = 200))
  # prob <- p200 * dmvnorm(gridData, mean = c(0, 0), sigma = sigmaAbn1) + p2X0 * dmvnorm(gridData, mean = c(muX, 0), sigma = sigmaAbn3) + 
  #   p20Y * dmvnorm(gridData, mean = c(0, muY), sigma = sigmaAbn3) + p2XY * dmvnorm(gridData, mean = c(muX, muY), sigma = sigmaAbn2)
  # sampData <- cbind(gridData, z = prob)
  # ggplot(sampData, aes(x = x, y = y, z = z, color = ..level..)) + 
  #   geom_contour() 
  
  dim(zk1) <- c(1, 2, K1, 1)
  dim(zk2) <- c(1, 2, K2, 1)
  
  simuData <- DfToRJafrocDataset(zk1, zk2, dataType = "ROC")
  desiredNumBins <- 5
  simuDataB <- DfBinDataset(simuData, desiredNumBins = desiredNumBins)
  
  ret <- FitCorCbmRoc(simuDataB, c(1, 1), c(1, 2), errBar = c(1, 4))
  c(ret$muX, ret$alphaX, ret$muY, ret$alphaY, ret$rhoNor, ret$rhoAbn2, ret$aucX, ret$aucY)
}
stopCluster(cl)

parmEstimates <- colMeans(ret)
parmNames <- c("muX", "alphaX", "muY", "alphaY", "rhoNor", "rhoAbn2", "aucX", "aucY")
names(parmEstimates) <- parmNames

aucX <- (1 - alphaX) * 0.5 + alphaX * pnorm(muX / sqrt(2))
aucY <- (1 - alphaY) * 0.5 + alphaY * pnorm(muY / sqrt(2))
parmTrue <- c(muX, alphaX, muY, alphaY, rhoNor, rhoAbn2, aucX, aucY)
parmDiff <- parmEstimates - parmTrue
print(parmDiff)

for (i in 1:6){
  tRet <- t.test(ret[ , i], mu = parmTrue[i], alternative = "t")
  cat("p-value for", parmNames[i], "is:", tRet$p.value, "\n")
}

i <- 7
tRet <- t.test(ret[ , i], mu = aucX, alternative = "t")
cat("p-value for", parmNames[i], "is:", tRet$p.value, "\n")

i <- 8
tRet <- t.test(ret[ , i], mu = aucY, alternative = "t")
cat("p-value for", parmNames[i], "is:", tRet$p.value, "\n")



