rm(list = ls()) # mainPlotAfroc.R

library(ggplot2)
library(grid)
library(RJafroc)

seed <- 1;set.seed(seed)
K1 <- 4;K2 <- 4
maxLL <- 2;Lk2 <- ceiling(runif(K2) * maxLL) 
mu <- 2;lambda <- 1;nu <- 1 ;zeta1 <- -1
frocDataRaw <- SimulateFrocDataset(mu = mu, lambda = lambda, nu = nu, I = 1, J = 1,
                             K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1)

maxNL <- dim(frocDataRaw$NL)[4]
FP <- apply(frocDataRaw$NL, 3, max)
FP <- FP[1:K1]
LL <- frocDataRaw$LL[(frocDataRaw$lesionWeight != -Inf)]
ratings <- c(FP,	LL)
truth <- c(rep(0, length(FP)), rep(1, length(LL)))
weight <- rep(1, length(ratings))
pos <- rep(0, length(ratings))
pos[1] <- 0.1
ret <- data.frame(ratings = ratings, truth = truth, weight = weight, pos = pos)
ret$ratings[ret$ratings == -Inf] <- -0.9
z_text <- data.frame(x = 3.5, y = -0.09, label = "z")

frocDataRaw$lesionWeight[3, ] <- c(0.6, 0.4)
frocDataRaw$lesionWeight[4, ] <- c(0.4, 0.6)
weight <- c(rep(1, length(FP)), frocDataRaw$lesionWeight[frocDataRaw$lesionWeight != -Inf])
ret <- data.frame(ratings = ratings, truth = truth, weight = weight, pos = pos)
ret$ratings[ret$ratings == -Inf] <- -0.9
z_text <- data.frame(x = 3.5, y = -0.09, label = "z")

fpScores <- sort(unique(ratings[truth == 0]), decreasing = TRUE)
llScores <- sort(unique(ratings[truth == 1]), decreasing = TRUE)
scores <- sort(unique(ratings), decreasing = TRUE)
fpf <- cumsum(as.vector(table(ratings[truth == 0])))/sum(truth == 0)
llf <- rep(NA, length(llScores))
for (z in 1:length(llScores)){
  llf[z] <- sum((ratings == llScores[z]) & (truth == 1))
}
llf <- cumsum(llf)/sum(truth == 1)
FPF <- 0:length(scores)
LLF <- FPF

numFP <- 1
numLL <- 1
for (k in 1:length(scores)) {
  if (!is.na(fpScores[numFP]) && fpScores[numFP] >= scores[k]) {
    FPF[k + 1] <- fpf[numFP]
    numFP <- numFP + 1
  } else {
    FPF[k + 1] <- FPF[k]
  }
  
  if (!is.na(llScores[numLL]) && llScores[numLL] >= scores[k]) {
    LLF[k + 1] <- llf[numLL]
    numLL <- numLL + 1
  } else {
    LLF[k + 1] <- LLF[k]
  }
}
FPF[k + 1] <- 1
LLF[k + 1] <- 1

AFROC_Points <- data.frame(FPF = FPF, LLF = LLF)
afrocPlot <- ggplot(AFROC_Points, aes(x = FPF, y = LLF) ) + 
  geom_line(size = 2) + geom_point(size = 5) + 
annotation_custom(grob = textGrob(bquote(italic("A")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.03, xmax = 0.03, 
                    ymin = 1/6 + 0.03, ymax = 1/6 + 0.03) +
  annotation_custom(grob = textGrob(bquote(italic("B")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.03, xmax = 0.03, 
                    ymin = 2/6 + 0.03, ymax = 2/6 + 0.03) +
  annotation_custom(grob = textGrob(bquote(italic("C")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.03, xmax = 0.03, 
                    ymin = 3/6 + 0.03, ymax = 3/6 + 0.03) +
  annotation_custom(grob = textGrob(bquote(italic("D")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.03, xmax = 0.03, 
                    ymin = 4/6 + 0.03, ymax = 4/6 + 0.03) +
  annotation_custom(grob = textGrob(bquote(italic("E")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.28, xmax = 0.28, 
                    ymin = 4/6 + 0.03, ymax = 4/6 + 0.03) +
  annotation_custom(grob = textGrob(bquote(italic("F")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.53, xmax = 0.53, 
                    ymin = 4/6 + 0.03, ymax = 4/6 + 0.03) +
  annotation_custom(grob = textGrob(bquote(italic("G")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.53, xmax = 0.53, 
                    ymin = 5/6 + 0.03, ymax = 5/6 + 0.03) +
  annotation_custom(grob = textGrob(bquote(italic("H")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.78, xmax = 0.78, 
                    ymin = 5/6 - 0.03, ymax = 5/6 - 0.03) +
  annotation_custom(grob = textGrob(bquote(italic("O")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = -0.03, xmax = -0.03, 
                    ymin = -0.03, ymax = -0.03) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) + scale_y_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0,1), x = c(0,1))
afrocPlot <- ggplotGrob(afrocPlot)
afrocPlot$layout$clip[afrocPlot$layout$name=="panel"] <- "off"
grid.draw(afrocPlot)
