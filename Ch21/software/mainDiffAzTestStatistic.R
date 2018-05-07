rm(list = ls()) # MainDiffAzTestStatistic.R

areaX  <-   0.8784
areaY  <-   0.8870
sigmaX  <-  0.0392     
sigmaY  <-  0.0376
rhoXY  <-   0.2569
varDiff <- sigmaX^2+sigmaY^2-2*rhoXY*sigmaX*sigmaY
testStat <- (areaX - areaY)/sqrt(varDiff) 
#testStat <- -1.644854 
cat("test statistic is :", testStat, "\n")
pValue2 <- pnorm(-abs(testStat))+(1-pnorm(abs(testStat)))
cat("two-tailed pValue, assuming AH: diff = 0, is :", pValue2, "\n")
pValue1 <- pnorm(testStat)
cat("one-tailed pValue is, assuming AH: diff < 0, is :", pValue1, "\n")
pValue1a <- 1 - pnorm(testStat)
cat("one-tailed pValue is, assuming AH: diff > 0, is :", pValue1a, "\n")
