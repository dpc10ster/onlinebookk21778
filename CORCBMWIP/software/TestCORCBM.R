rm(list = ls())
library(RJafroc)
library(ggplot2)


library(mvtnorm)
library(binom)

seed <- 123;set.seed(seed)
muX <- 1.5;muY <- 3
alphaX <- 0.4;alphaY <- 0.7
rhoNor <- 0.3;rhoAbn1 <- rhoNor;rhoAbn2 <- 0.8
rhoAbn3 <- mean(c(rhoAbn1, rhoAbn2))
sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
sigmaAbn1 <- rbind(c(1, rhoAbn1), c(rhoAbn1, 1))
sigmaAbn2 <- rbind(c(1, rhoAbn2), c(rhoAbn2, 1))
sigmaAbn3 <- rbind(c(1, rhoAbn3), c(rhoAbn3, 1))
# 50/50; 100/100; 1000/1000; 5000/5000
K1 <- 51;K2 <- 50
p200 <- (1 - alphaX) * (1 - alphaY)
p2X0 <- alphaX * (1 - alphaY)
p20Y <- (1 - alphaX) * alphaY
p2XY <- alphaX * alphaY
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

dim(zk1) <- c(1, 2, K1)
dim(zk2) <- c(1, 2, K2)

simuData <- DfToRJafrocDataset(zk1, zk2)

ret <- FitCorCbmRoc(simuData)