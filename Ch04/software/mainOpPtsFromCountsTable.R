# mainOpPtsFromCountsTable.R
rm( list = ls())
library(ggplot2)
source("plotROC.R")

options(digits = 6)
Ktr = array(dim = c(2,5))
# Table 4.2.1
Ktr[1,]  <- c(30,19,8,2,1) 
Ktr[2,]  <- c(5,6,5,12,22)

R <- length(Ktr[1,]) -1

FPF <- array(0,dim = R)
TPF <- array(0,dim = R)

for (r in (R+1):2) {
  FPF[(R+2)-r] <- 
    sum(Ktr[1, r:(R+1)])/sum(Ktr[1,])
  TPF[(R+2)-r] <- 
    sum(Ktr[2, r:(R+1)])/sum(Ktr[2,])    
}

cat("FPF =", "\n")
cat(FPF, "\n")
cat("TPF = ", "\n")
cat(TPF, "\n")

mu <- qnorm(.5)+qnorm(.9)
sigma <- 1
Az <- pnorm(mu/sqrt(2))

zeta <- seq(- mu - 2,+ mu + 2,0.01)
fpf <- array(dim = length(zeta))
tpf <- array(dim = length(zeta))
for (i in 1:length(zeta)) {
  fpf[i] <- pnorm(-zeta[i])
  tpf[i] <- pnorm((mu -zeta[i])/sigma) 
}
curveData <- data.frame(FPF = c(1, fpf, 0), TPF = c(1, tpf, 0))
pointsData <- data.frame(FPF = FPF, TPF = TPF)
rocPlot1 <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = curveData, size = 2) + 
  geom_point(data = pointsData, size = 5) + 
  xlab("FPF")+ ylab("TPF" ) + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
print(rocPlot1)

cat("uppermost point based estimate of mu = ", mu, "\n")
cat("corresponding estimate of Az = ", Az, "\n")

mu <- 2.17;sigma <- 1.65
Az <- pnorm(mu/sqrt(1+sigma^2))
zeta <- seq(- mu - 2,+ mu + 2,0.01)
fpf <- array(dim = length(zeta))
tpf <- array(dim = length(zeta))
for (i in 1:length(zeta)) {
  fpf[i] <- pnorm(-zeta[i])
  tpf[i] <- pnorm((mu -zeta[i])/sigma) 
}

zeta <- seq(- mu - 2,+ mu + 2,0.01)
fpf <- array(dim = length(zeta))
tpf <- array(dim = length(zeta))
for (i in 1:length(zeta)) {
  fpf[i] <- pnorm(-zeta[i])
  tpf[i] <- pnorm((mu -zeta[i])/sigma) 
}
curveData <- data.frame(FPF = c(1, fpf, 0), TPF = c(1, tpf, 0))
pointsData <- data.frame(FPF = FPF, TPF = TPF)
rocPlot2 <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = curveData, size = 2) + 
  geom_point(data = pointsData, size = 5) + 
  xlab("FPF") + ylab("TPF" ) + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
print(rocPlot2)

cat("binormal estimate of Az = ", Az, "\n")
