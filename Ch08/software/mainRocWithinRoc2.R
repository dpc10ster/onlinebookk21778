rm(list = ls()) # mainRocWithinRoc.R
require("ggplot2");require("grid")
source("Wilcoxon.R")

seed <- 1;set.seed(seed)
muNH <- 1.5;muAH <- 2.1;sigma <- 1.3
K1 <- 50;K2 <- 52#;K1 <- K1*2;K2 <- K2*2

# cheat to find the population mean and std. dev.
AUC <- array(dim = 10000)
for (i in 1:length(AUC)) {
  zk1 <- rnorm(K1);zk2 <- rnorm(K2, mean = muNH, sd = sigma)  
  AUC[i] <- Wilcoxon(zk1, zk2)
}
sigmaAUC <- sqrt(var(AUC));meanAUC <- mean(AUC)

T <- 2000
mu <- c(muNH,muAH)
alphaArr <- seq(0.05, 0.95, length.out = 10)
EmpAlpha <- array(dim = length(alphaArr))
EmpPower <- array(dim = length(alphaArr))
for (a in 1:length(alphaArr)) {
  alpha <- alphaArr[a]
  reject <- array(0, dim = c(2, T))
  for (h in 1:2) {  
    for (t in 1:length(reject[h,])) {  
      zk1 <- rnorm(K1);zk2 <- rnorm(K2, mean = mu[h], sd = sigma)  
      AUC <- Wilcoxon(zk1, zk2)  
      obsvdZ <- (AUC - meanAUC)/sigmaAUC
      p <- 2*pnorm(-abs(obsvdZ)) # p value for individual t
      if (p < alpha) reject[h,t] = 1 
    }
  }
  EmpAlpha[a] <- sum(reject[1,])/length(reject[1,])
  EmpPower[a] <- sum(reject[2,])/length(reject[2,])
}
# EmpAlpha <- c(0,EmpAlpha,1)
# EmpPower <- c(0,EmpPower,1)

pointData <- data.frame(EmpAlpha = EmpAlpha, EmpPower = EmpPower)
zetas <- seq(-5, 5, by = 0.01)
muRoc <- 1.8
curveData <- data.frame(EmpAlpha = pnorm(-zetas), 
                        EmpPower = pnorm(muRoc - zetas))
p <- ggplot(mapping = aes(x = EmpAlpha, y = EmpPower)) + 
  geom_point(data = pointData, size = 5) + 
  geom_line(data = curveData, size = 2) + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) +
  annotation_custom(grob = textGrob(bquote(italic("O")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = -0.03, xmax = -0.03, # adjust the position of "O"
                    ymin = -0.03, ymax = -0.03) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0,1), x = c(0,1))  
  
p <- ggplotGrob(p)
p$layout$clip[p$layout$name=="panel"] <- "off"
grid.draw(p)
