rm(list = ls()) # mainProofPlot.R
library(ggplot2)
library(grid)
library(pBrackets)
ratings <- c(2.012014287,	2.303861492,	-2.67621502,	-0.731465008,	
             2.012014287,	2.923751322, 0.778163587, -1.080014238,	-2.67621502)
truth <- c(rep(0, 4), rep(1, 5))
weight <- c(rep(5, 4), 3, 5, 2, 2, 3) / 5
pos <- rep(0, 9)
pos[5] <- 0.1
ret <- data.frame(ratings = ratings, truth = truth, weight = weight, pos = pos)

fpScores <- sort(unique(ratings[truth == 0]), decreasing = TRUE)
llScores <- sort(unique(ratings[truth == 1]), decreasing = TRUE)
scores <- sort(unique(ratings), decreasing = TRUE)
fpf <- cumsum(as.vector(table(ratings[truth == 0])))/sum(truth == 0)
wllf <- rep(NA, length(llScores))
for (z in 1:length(llScores)){
  wllf[z] <- sum(weight[(ratings == llScores[z]) & (truth == 1)])
}
wllf <- cumsum(wllf)/sum(weight[5:9])
FPF <- 0:length(scores)
wLLF <- FPF

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
    wLLF[k + 1] <- wllf[numLL]
    numLL <- numLL + 1
  } else {
    wLLF[k + 1] <- wLLF[k]
  }
}
FPF[k + 1] <- 1
wLLF[k + 1] <- 1

wAFROC_Points <- data.frame(FPF = FPF, wLLF = wLLF)

zetaZero_x <- wAFROC_Points$FPF[length(wAFROC_Points$FPF)]
zetaZero_y <- wAFROC_Points$wLLF[length(wAFROC_Points$wLLF)]
Zero <- data.frame(x = zetaZero_x, y = zetaZero_y, label = "0")
zetaOne_x <- wAFROC_Points$FPF[length(wAFROC_Points$FPF) - 1]
zetaOne_y <- wAFROC_Points$wLLF[length(wAFROC_Points$wLLF) - 1]
One <- data.frame(x = zetaOne_x, y = zetaOne_y, label = "1")
zeta1_x <- wAFROC_Points$FPF[3]
zeta1_y <- wAFROC_Points$wLLF[3]
zeta1 <- data.frame(x = c(zeta1_x, zeta1_x), y = c(0, zeta1_y))
i1 <- data.frame(x = zeta1_x, y = zeta1_y, label = "i + 1")
zeta2_x <- wAFROC_Points$FPF[4]
zeta2_y <- wAFROC_Points$wLLF[4]
zeta2 <- data.frame(x = c(zeta2_x, zeta2_x), y = c(0, zeta2_y))
i2 <- data.frame(x = zeta2_x, y = zeta2_y, label = "i")

shade1 <- data.frame(FPF = c(zeta2_x, zeta1_x, zeta1_x, zeta2_x), wLLF = c(zeta2_y, zeta1_y, 0, 0))
shade2 <- data.frame(FPF = c(1, 1, 0.75, 0.75), wLLF = c(1, 0, 0, 0.8))

p <- ggplot(wAFROC_Points, aes(x = FPF, y = wLLF) ) + 
  geom_polygon(data = shade1, fill = 'grey') + 
  geom_polygon(data = shade2, fill = 'grey') + 
  geom_line(size = 2) + geom_point(size = 5) + 
  geom_line(data = zeta1, aes(x = x, y = y), linetype = 2) + geom_line(data = zeta2, aes(x = x, y = y), linetype = 2) +
  geom_segment(aes(x = 0.75, y = 0, xend = 0.75, yend = 0.8), linetype = 2) +  
  geom_segment(aes(x = 1, y = 0.8, xend = 0.75, yend = 0.8), linetype = 2) +  
  geom_text(data = Zero, aes(x = x, y = y, label = label), vjust = 2, hjust = 1.1, size = 9) + 
  geom_text(data = One, aes(x = x, y = y, label = label), vjust = 1.5, hjust = -0.5, size = 10) + 
  geom_text(data = i1, aes(x = x, y = y, label = label), vjust = -0.5, hjust = 1, size = 10) + 
  geom_text(data = i2, aes(x = x, y = y, label = label), vjust = -0.5, hjust = 3, size = 10) + 
  annotation_custom(grob = textGrob(bquote(italic(wLLF["1"])),
                                    gp = gpar(fontsize = 22)), 
                    xmin = 0.87, xmax = 0.87, 
                    ymin = 0.47, ymax = 0.47) +
  geom_segment(aes(x = 0.85, y = 0.5, xend = 0.75, yend = 0.6), arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  annotation_custom(grob = textGrob(bquote(italic(wLLF[i + "1"])),
                                    gp = gpar(fontsize = 22)), 
                    xmin = zeta1_x - 0.15, xmax = zeta1_x - 0.05, 
                    ymin = 0.25, ymax = 0.25) +
  geom_segment(aes(x = zeta1_x - 0.1, y = 0.22, xend = zeta1_x, yend = 0.15), arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  annotation_custom(grob = grobTree(textGrob(bquote(italic(wLLF[i])),
                                             gp = gpar(fontsize = 22))), 
                    xmin = zeta2_x + 0.05, xmax = zeta2_x + 0.15, 
                    ymin = 0.35, ymax = 0.35) +
  geom_segment(aes(x = zeta2_x + 0.1, y = 0.32, xend = zeta2_x, yend = 0.25), arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  annotation_custom(grob = textGrob(bquote(italic(FPF[i]-FPF[i + "1"])),
                                    gp = gpar(fontsize = 18)), 
                    xmin = 0.375, xmax = 0.375, 
                    ymin = 0.1, ymax = 0.1) +
  geom_segment(aes(x = 0.375, y = 0.09, xend = 0.3, yend = 0), arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  annotation_custom(grob = textGrob(bquote(italic(1-FPF["1"])),
                                    gp = gpar(fontsize = 18)), 
                    xmin = 0.875, xmax = 0.875, 
                    ymin = 0.7, ymax = 0.7) +
  geom_segment(aes(x = 0.875, y = 0.72, xend = 0.9, yend = 0.8), arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text = element_text(size = 15), axis.title = element_text(size = 20)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0,1), x = c(0,1)) +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))

        print(p)
        