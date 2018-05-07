# mainIsFrocGood.R
rm(list = ls())
library(RJafroc)
library(ggplot2)

logseq <- function( d1, d2, n) {
  logf <- log(d2/d1)/(n-1)
  return (exp(seq(log(d1), log(d2), logf)))
}

Lmax <- 1;K2 <- 700;Lk2 <- floor(runif(K2, 1, Lmax + 1))
nLesPerCase <- unique(Lk2)
lesionDist <- array(dim = c(length(nLesPerCase), 2))
for (i in nLesPerCase) lesionDist[i, ] <- c(i, sum(Lk2 == i)/K2)

## PART I
lambda <- 1;nu <- 1
cat("Vary mu only, lambda and nu equal to 1\n")
muArr <- logseq(0.001, 10, 9)
aucRoc <- rep(NA, length(muArr))
aucAfroc <- aucRoc;aucFroc <- aucRoc
for (i in 1:length(muArr)) { # vary mu loop
  mu <- muArr[i]
  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu,
    lesionDistribution = lesionDist
  )
  cat("mu = ", mu,", lambda = ", lambda,
      ", nu = ", nu, ", aucFroc = ", ret1$aucFROC,
      ", aucRoc = ", ret1$aucROC,", aucAfroc = ", ret1$aucAFROC,"\n")
  aucRoc[i] <- ret1$aucROC
  aucAfroc[i] <- ret1$aucAFROC
  aucFroc[i] <- ret1$aucFROC
}
cat("approx slope AFROC vs ROC =", 
    (aucAfroc[9]-aucAfroc[1])/(aucRoc[9]-aucRoc[1]),"\n")
aucRocFroc <- data.frame(aucRoc = aucRoc, aucFroc = aucFroc)
aucRocAfroc <- data.frame(aucRoc = aucRoc, aucAfroc = aucAfroc)

muArr2 <- seq(1, 10, by = 0.1)
aucRoc2 <- rep(NA, length(muArr))
aucAfroc2 <- aucRoc2;aucFroc2 <- aucRoc2
for (i in 1:length(muArr2)) { # vary mu loop
  mu <- muArr2[i]
  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu,
    lesionDistribution = lesionDist
  )
  aucRoc2[i] <- ret1$aucROC
  aucAfroc2[i] <- ret1$aucAFROC
  aucFroc2[i] <- ret1$aucFROC
}
plotCurveRocFroc <- data.frame(
  aucRoc = c(aucRoc, aucRoc2), aucFroc = c(aucFroc, aucFroc2))
plotRocFroc <- ggplot(
  data = aucRocFroc, mapping = aes(x = aucRoc, y = aucFroc)) + 
  geom_point(size = 5) + 
  geom_line(data = plotCurveRocFroc, size = 2) + 
  xlab("ROC-AUC") + ylab("FROC-AUC")  +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(plotRocFroc)

plotCurveRocAfroc <- data.frame(
  aucRoc = c(aucRoc, aucRoc2), aucAfroc = c(aucAfroc, aucAfroc2))
plotRocAfroc <- ggplot(
  data = aucRocAfroc, mapping = aes(x = aucRoc, y = aucAfroc)) + 
  geom_point(size = 5) + 
  geom_line(data = plotCurveRocAfroc, size = 2) + 
  xlab("ROC-AUC") + ylab("AFROC-AUC")   +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(plotRocAfroc)

## PART II
mu <- 2;nu <- 1
cat("\nVary lambda only, mu = 2 and nu = 1\n")
lambdaArr <- logseq(0.2, 8, 9)
aucRoc <- rep(NA, length(lambdaArr))
aucAfroc <- aucRoc;aucFroc <- aucRoc
for (i in 1:length(lambdaArr)) {
  lambda <- lambdaArr[i]
  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, 
    lesionDistribution = lesionDist
  ) 
  cat("mu = ", mu,", lambda = ", lambda,
      ", nu = ", nu, ", aucFroc = ", ret1$aucFROC,
      ", aucRoc = ", ret1$aucROC,", aucAfroc = ", ret1$aucAFROC,"\n")
  aucRoc[i] <- ret1$aucROC
  aucAfroc[i] <- ret1$aucAFROC
  aucFroc[i] <- ret1$aucFROC
}
cat("approx slope AFROC vs ROC =", 
    (aucAfroc[9]-aucAfroc[1])/(aucRoc[9]-aucRoc[1]),"\n")
aucRocFroc <- data.frame(aucRoc = aucRoc, aucFroc = aucFroc)
aucRocAfroc <- data.frame(aucRoc = aucRoc, aucAfroc = aucAfroc)

lambdaArr2 <- seq(0.2, 8, by = 0.1)
aucRoc2 <- rep(NA, length(lambdaArr2))
aucAfroc2 <- aucRoc2;aucFroc2 <- aucRoc2
for (i in 1:length(lambdaArr2)) { # vary mu loop
  lambda <- lambdaArr2[i]
  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu,
    lesionDistribution = lesionDist
  )
  aucRoc2[i] <- ret1$aucROC
  aucAfroc2[i] <- ret1$aucAFROC
  aucFroc2[i] <- ret1$aucFROC
}
plotCurveRocFroc <- data.frame(
  aucRoc = c(aucRoc, aucRoc2), aucFroc = c(aucFroc, aucFroc2))
plotRocFroc <- ggplot(
  data = aucRocFroc, mapping = aes(x = aucRoc, y = aucFroc)) + 
  geom_point(size = 5) + 
  geom_line(data = plotCurveRocFroc, size = 2) + 
  xlab("ROC-AUC") + ylab("FROC-AUC") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(plotRocFroc)

plotCurveRocAfroc <- data.frame(
  aucRoc = c(aucRoc, aucRoc2), aucAfroc = c(aucAfroc, aucAfroc2))
plotRocAfroc <- ggplot(
  data = aucRocAfroc, mapping = aes(x = aucRoc, y = aucAfroc)) + 
  geom_point(size = 5) + 
  geom_line(data = plotCurveRocAfroc, size = 2) + 
  xlab("ROC-AUC") + ylab("AFROC-AUC") +
  theme(axis.title.y = 
          element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(plotRocAfroc)

## PART III
mu <- 2;lambda <- 1
cat("\nVary nu only, mu = 2 and lambda = 1\n")
nuArr <- logseq(0.2, 8, 9)
aucRoc <- rep(NA, length(nuArr));aucAfroc <- aucRoc;aucFroc <- aucRoc
for (i in 1:length(nuArr)) {
  nu <- nuArr[i]
  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, 
    lesionDistribution = lesionDist,
    llfRange = c(0,1)
  ) 
  cat("mu = ", mu,", lambda = ", lambda,
      ", nu = ", nu, ", aucFroc = ", ret1$aucFROC,
      ", aucRoc = ", ret1$aucROC,", aucAfroc = ", ret1$aucAFROC,"\n")
  aucRoc[i] <- ret1$aucROC
  aucAfroc[i] <- ret1$aucAFROC
  aucFroc[i] <- ret1$aucFROC
}
cat("approx slope AFROC vs ROC =", 
    (aucAfroc[9]-aucAfroc[1])/(aucRoc[9]-aucRoc[1]),"\n")
aucRocFroc <- data.frame(aucRoc = aucRoc, aucFroc = aucFroc)
aucRocAfroc <- data.frame(aucRoc = aucRoc, aucAfroc = aucAfroc)

nuArr2 <- seq(0.2, 8, by = 0.1)
aucRoc2 <- rep(NA, length(nuArr2))
aucAfroc2 <- aucRoc2;aucFroc2 <- aucRoc2
for (i in 1:length(nuArr2)) { # vary mu loop
  nu <- nuArr[2]
  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu,
    lesionDistribution = lesionDist
  )
  aucRoc2[i] <- ret1$aucROC
  aucAfroc2[i] <- ret1$aucAFROC
  aucFroc2[i] <- ret1$aucFROC
}
plotCurveRocFroc <- data.frame(
  aucRoc = c(aucRoc, aucRoc2), aucFroc = c(aucFroc, aucFroc2))
plotRocFroc <- ggplot(
  data = aucRocFroc, mapping = aes(x = aucRoc, y = aucFroc)) + 
  geom_point(size = 5) + 
  geom_line(data = plotCurveRocFroc, size = 2) + 
  xlab("ROC-AUC") + ylab("FROC-AUC") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(plotRocFroc)
plotCurveRocAfroc <- data.frame(
  aucRoc = c(aucRoc, aucRoc2), aucAfroc = c(aucAfroc, aucAfroc2))
plotRocAfroc <- ggplot(
  data = aucRocAfroc, mapping = aes(x = aucRoc, y = aucAfroc)) + 
  geom_point(size = 5) + 
  geom_line(data = plotCurveRocAfroc, size = 2) + 
  xlab("ROC-AUC") + ylab("AFROC-AUC") +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(plotRocAfroc)