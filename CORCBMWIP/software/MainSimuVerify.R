rm(list = ls())
library(mvtnorm)
library(RJafroc)
library(ggplot2)
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
K1 <- 50;K2 <- 50
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

simuData <- Df2RJafrocDataset(zk1, zk2)
simuDataB <- DfBinDataset(simuData, desiredNumBins = 5)

ret <- FitCorCbm(simuDataB)

fittedPlot <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = ret$fittedPlot$plot_env$plotCBM, mapping = aes(linetype = Condition), size = 1) + 
  geom_point(data = ret$fittedPlot$plot_env$plotOpPnts, mapping = aes(shape = Condition), size = 3) + 
  theme(legend.title=element_blank(), legend.position = c(1, 0), legend.direction = "horizontal", legend.justification = c(1, 0), legend.key.size = unit(1, "cm")) + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

ciIndx <- c(1, 4)
FPF <- ret$fittedPlot$plot_env$FPFX[ciIndx]
TPF <- ret$fittedPlot$plot_env$TPFX[ciIndx]
ciX <- binom.confint(x = FPF * K1, n = K1, method = "exact")
ciY <- binom.confint(x = TPF * K2, n = K2, method = "exact")
ciXUpper <- ciX$upper
ciXLower <- ciX$lower
ciYUpper <- ciY$upper
ciYLower <- ciY$lower
for (p in 1:length(FPF)){
  ciX <- data.frame(FPF = c(ciXUpper[p], ciXLower[p]), TPF = c(TPF[p], TPF[p]))
  ciY <- data.frame(FPF = c(FPF[p], FPF[p]), TPF = c(ciYUpper[p], ciYLower[p]))
  fittedPlot <- fittedPlot + geom_line(data = ciY, aes(x = FPF, y = TPF), color = "black") + 
    geom_line(data = ciX, aes(x = FPF, y = TPF), color = "black")
  barRgt <- data.frame(FPF = c(ciXUpper[p], ciXUpper[p]), TPF = c(TPF[p] - 0.01, TPF[p] + 0.01))
  barLft <- data.frame(FPF = c(ciXLower[p], ciXLower[p]), TPF = c(TPF[p] - 0.01, TPF[p] + 0.01))
  barUp <- data.frame(FPF = c(FPF[p] - 0.01, FPF[p] + 0.01), TPF = c(ciYUpper[p], ciYUpper[p]))
  barBtm <- data.frame(FPF = c(FPF[p] - 0.01, FPF[p] + 0.01), TPF = c(ciYLower[p], ciYLower[p]))
  fittedPlot <- fittedPlot + geom_line(data = barRgt, aes(x = FPF, y = TPF), color = "black") + 
    geom_line(data = barLft, aes(x = FPF, y = TPF), color = "black") + 
    geom_line(data = barUp, aes(x = FPF, y = TPF), color = "black") + 
    geom_line(data = barBtm, aes(x = FPF, y = TPF), color = "black")
}

FPF <- ret$fittedPlot$plot_env$FPFY[ciIndx]
TPF <- ret$fittedPlot$plot_env$TPFY[ciIndx]
ciX <- binom.confint(x = FPF * K1, n = K1, method = "exact")
ciY <- binom.confint(x = TPF * K2, n = K2, method = "exact")
ciXUpper <- ciX$upper
ciXLower <- ciX$lower
ciYUpper <- ciY$upper
ciYLower <- ciY$lower
for (p in 1:length(FPF)){
  ciX <- data.frame(FPF = c(ciXUpper[p], ciXLower[p]), TPF = c(TPF[p], TPF[p]))
  ciY <- data.frame(FPF = c(FPF[p], FPF[p]), TPF = c(ciYUpper[p], ciYLower[p]))
  fittedPlot <- fittedPlot + geom_line(data = ciY, aes(x = FPF, y = TPF), color = "black") + 
    geom_line(data = ciX, aes(x = FPF, y = TPF), color = "black")
  barRgt <- data.frame(FPF = c(ciXUpper[p], ciXUpper[p]), TPF = c(TPF[p] - 0.01, TPF[p] + 0.01))
  barLft <- data.frame(FPF = c(ciXLower[p], ciXLower[p]), TPF = c(TPF[p] - 0.01, TPF[p] + 0.01))
  barUp <- data.frame(FPF = c(FPF[p] - 0.01, FPF[p] + 0.01), TPF = c(ciYUpper[p], ciYUpper[p]))
  barBtm <- data.frame(FPF = c(FPF[p] - 0.01, FPF[p] + 0.01), TPF = c(ciYLower[p], ciYLower[p]))
  fittedPlot <- fittedPlot + geom_line(data = barRgt, aes(x = FPF, y = TPF), color = "black") + 
    geom_line(data = barLft, aes(x = FPF, y = TPF), color = "black") + 
    geom_line(data = barUp, aes(x = FPF, y = TPF), color = "black") + 
    geom_line(data = barBtm, aes(x = FPF, y = TPF), color = "black")
}
print(fittedPlot)

ciHalfWidth <- ret$stdErr * qnorm(0.975)
cat("Original muX is:", muX, "estimate is", ret$muX, ", 95% CI is", ret$muX - ciHalfWidth[1], ",", ret$muX + ciHalfWidth[1], 
    ". Included?", ifelse((ret$muX - ciHalfWidth[1] - muX) * (ret$muX + ciHalfWidth[1] - muX) <= 0, yes = "Yes", no = "No"), "\n")

cat("Original alphaX is:", alphaX, "estimate is", ret$alphaX, ", 95% CI is", ret$alphaX - ciHalfWidth[2], ",", ret$alphaX + ciHalfWidth[2], 
    ". Included?", ifelse((ret$alphaX - ciHalfWidth[2] - alphaX) * (ret$alphaX + ciHalfWidth[2] - alphaX) <= 0, yes = "Yes", no = "No"), "\n")

cat("Original muY is:", muY, "estimate is", ret$muY, ", 95% CI is", ret$muY - ciHalfWidth[3], ",", ret$muY + ciHalfWidth[3], 
    ". Included?", ifelse((ret$muY - ciHalfWidth[3] - muY) * (ret$muY + ciHalfWidth[3] - muY) <= 0, yes = "Yes", no = "No"), "\n")

cat("Original alphaY is:", alphaY, "estimate is", ret$alphaY, ", 95% CI is", ret$alphaY - ciHalfWidth[4], ",", ret$alphaY + ciHalfWidth[4], 
    ". Included?", ifelse((ret$alphaY - ciHalfWidth[4] - alphaY) * (ret$alphaY + ciHalfWidth[4] - alphaY) <= 0, yes = "Yes", no = "No"), "\n")

cat("Original rhoNor is:", rhoNor, "estimate is", ret$rhoNor, ", 95% CI is", ret$rhoNor - ciHalfWidth[5], ",", ret$rhoNor + ciHalfWidth[5], 
    ". Included?", ifelse((ret$rhoNor - ciHalfWidth[5] - rhoNor) * (ret$rhoNor + ciHalfWidth[5] - rhoNor) <= 0, yes = "Yes", no = "No"), "\n")

cat("Original rhoAbn2 is:", rhoAbn2, "estimate is", ret$rhoAbn2, ", 95% CI is", ret$rhoAbn2 - ciHalfWidth[6], ",", ret$rhoAbn2 + ciHalfWidth[6], 
    ". Included?", ifelse((ret$rhoAbn2 - ciHalfWidth[6] - rhoAbn2) * (ret$rhoAbn2 + ciHalfWidth[6] - rhoAbn2) <= 0, yes = "Yes", no = "No"), "\n")

cat(ret$muX, ret$muY, ret$alphaX, ret$alphaY, ret$rhoNor, ret$rhoAbn2, ret$aucX, ret$aucY, "\n")
cat(ret$muX - ciHalfWidth[1], ",", ret$muX + ciHalfWidth[1], ret$muY - ciHalfWidth[3], ",", ret$muY + ciHalfWidth[3], ret$alphaX - ciHalfWidth[2], ",", ret$alphaX + ciHalfWidth[2], 
    ret$alphaY - ciHalfWidth[4], ",", ret$alphaY + ciHalfWidth[4], ret$rhoNor - ciHalfWidth[5], ",", ret$rhoNor + ciHalfWidth[5], ret$rhoAbn2 - ciHalfWidth[6], ",", ret$rhoAbn2 + ciHalfWidth[6],
    ret$aucX - ciHalfWidth[5], ",", ret$aucX + ciHalfWidth[5], ret$aucY - ciHalfWidth[6], ",", ret$aucY + ciHalfWidth[6])

# > source('~/xzwps/CorCBM/software/MainSimuVerify.R')
# Loading required package: shiny
# Loading required package: xlsx
# Loading required package: rJava
# Loading required package: xlsxjars
# Original muX is: 1.5 estimate is 1.106443 , 95% CI is -0.4303567 , 2.643243 . Included? Yes 
# Original alphaX is: 0.4 estimate is 0.595083 , 95% CI is -0.1621231 , 1.352289 . Included? Yes 
# Original muY is: 3 estimate is 3.997848 , 95% CI is 1.537069 , 6.458627 . Included? Yes 
# Original alphaY is: 0.7 estimate is 0.7234784 , 95% CI is 0.5896488 , 0.857308 . Included? Yes 
# Original rhoNor is: 0.3 estimate is 0.3337149 , 95% CI is 0.04237823 , 0.6250515 . Included? Yes 
# Original rhoAbn2 is: 0.8 estimate is 0.4455214 , 95% CI is -1.254607 , 2.14565 . Included? Yes 
# 1.106443 3.997848 0.595083 0.7234784 0.3337149 0.4455214 0.6684099 0.860039 
# -0.4303567 , 2.643243 1.537069 , 6.458627 -0.1621231 , 1.352289 0.5896488 , 0.857308 0.04237823 , 0.6250515 -1.254607 , 2.14565 0.3770733 , 0.9597466 -0.8400897 , 2.560168
# > source('~/xzwps/CorCBM/software/MainSimuVerify.R')
# Original muX is: 1.5 estimate is 0.8516743 , 95% CI is -0.0002817741 , 1.70363 . Included? Yes 
# Original alphaX is: 0.4 estimate is 0.8154062 , 95% CI is 0.1182744 , 1.512538 . Included? Yes 
# Original muY is: 3 estimate is 3.222437 , 95% CI is 2.586523 , 3.858352 . Included? Yes 
# Original alphaY is: 0.7 estimate is 0.6846826 , 95% CI is 0.5845043 , 0.784861 . Included? Yes 
# Original rhoNor is: 0.3 estimate is 0.1875667 , 95% CI is -0.03079933 , 0.4059326 . Included? Yes 
# Original rhoAbn2 is: 0.8 estimate is 0.8035768 , 95% CI is 0.4132559 , 1.193898 . Included? Yes 
# 0.8516743 3.222437 0.8154062 0.6846826 0.1875667 0.8035768 0.6846795 0.8345735 
# -0.0002817741 , 1.70363 2.586523 , 3.858352 0.1182744 , 1.512538 0.5845043 , 0.784861 -0.03079933 , 0.4059326 0.4132559 , 1.193898 0.4663135 , 0.9030455 0.4442526 , 1.224894
# > source('~/xzwps/CorCBM/software/MainSimuVerify.R')
# Original muX is: 1.5 estimate is 1.884681 , 95% CI is 1.271864 , 2.497498 . Included? Yes 
# Original alphaX is: 0.4 estimate is 0.3346187 , 95% CI is 0.2462637 , 0.4229736 . Included? Yes 
# Original muY is: 3 estimate is 3.151775 , 95% CI is 2.852606 , 3.450944 . Included? Yes 
# Original alphaY is: 0.7 estimate is 0.6950135 , 95% CI is 0.6601503 , 0.7298767 . Included? Yes 
# Original rhoNor is: 0.3 estimate is 0.2765468 , 95% CI is 0.215264 , 0.3378296 . Included? Yes 
# Original rhoAbn2 is: 0.8 estimate is 0.908306 , 95% CI is 0.5822487 , 1.234363 . Included? Yes 
# 1.884681 3.151775 0.3346187 0.6950135 0.2765468 0.908306 0.636752 0.8385279 
# 1.271864 , 2.497498 2.852606 , 3.450944 0.2462637 , 0.4229736 0.6601503 , 0.7298767 0.215264 , 0.3378296 0.5822487 , 1.234363 0.5754692 , 0.6980348 0.5124707 , 1.164585
# > source('~/xzwps/CorCBM/software/MainSimuVerify.R')
# Original muX is: 1.5 estimate is 1.609413 , 95% CI is 1.371434 , 1.847392 . Included? Yes 
# Original alphaX is: 0.4 estimate is 0.3659602 , 95% CI is 0.3178065 , 0.4141139 . Included? Yes 
# Original muY is: 3 estimate is 3.01184 , 95% CI is 2.851064 , 3.172616 . Included? Yes 
# Original alphaY is: 0.7 estimate is 0.7028191 , 95% CI is 0.6848662 , 0.7207721 . Included? Yes 
# Original rhoNor is: 0.3 estimate is 0.2924337 , 95% CI is 0.2649134 , 0.3199541 . Included? Yes 
# Original rhoAbn2 is: 0.8 estimate is 0.6079936 , 95% CI is 0.3611145 , 0.8548727 . Included? Yes 
# 1.609413 3.01184 0.3659602 0.7028191 0.2924337 0.6079936 0.6363002 0.8397438 
# 1.371434 , 1.847392 2.851064 , 3.172616 0.3178065 , 0.4141139 0.6848662 , 0.7207721 0.2649134 , 0.3199541 0.3611145 , 0.8548727 0.6087798 , 0.6638206 0.5928647 , 1.086623