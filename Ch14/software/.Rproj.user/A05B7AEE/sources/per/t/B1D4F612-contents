rm(list = ls()) # mainLinearPlotAfroc.R

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

linearPlotAFROC <- ggplot(data = ret) + geom_point(aes(x = ratings, y = pos, color = factor(truth)), size = 8) + 
  scale_color_manual(values = c("green", "red")) + 
  scale_size_manual(values = c(2, 3, 4, 5) + 3) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), legend.position="none") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) +
  scale_y_continuous(name = "", breaks = NULL) + scale_x_continuous(name = "", breaks = NULL) + 
  geom_segment(aes(x = -1, y = -0.05, xend = 3.4, yend = -0.05), arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  geom_text(data = z_text, aes(x = x, y = y, label = label), vjust = -0.5, hjust = 1, size = 10) + 
  geom_segment(aes(x = 3.3, y = -0.4, xend = 3.3, yend = 0.4), linetype = 2) +
  geom_segment(aes(x = 3.3, y = 0.25, xend = 3.15, yend = 0.25), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(zeta),
                                    gp = gpar(fontsize = 30)), 
                    xmin = 3.2, xmax = 3.2, 
                    ymin = 0.35, ymax = 0.35) +
  geom_segment(aes(x = -0.9, y = -0.4, xend = -0.9, yend = 0.4), linetype = 2) +
  annotation_custom(grob = textGrob(bquote(z == -infinity),
                                    gp = gpar(fontsize = 25)), 
                    xmin = -0.9, xmax = -0.9, 
                    ymin = -0.42, ymax = -0.42) +
  geom_segment(aes(x = -0.7, y = 0.25, xend = -0.87, yend = 0.12), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(FP[11]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = -0.68, xmax = -0.68, 
                    ymin = 0.28, ymax = 0.28) +
  geom_segment(aes(x = -0.45, y = 0.15, xend = -0.33, yend = 0.03), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(FP[41]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = -0.42, xmax = -0.42, 
                    ymin = 0.18, ymax = 0.18) +
  geom_segment(aes(x = 0.33, y = 0.15, xend = 0.45, yend = 0.03), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(FP[21]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = 0.33, xmax = 0.33, 
                    ymin = 0.19, ymax = 0.19) +
  geom_segment(aes(x = 0.9, y = 0.15, xend = 0.76, yend = 0.03), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(FP[31]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = 0.9, xmax = 0.9, 
                    ymin = 0.19, ymax = 0.19) +
  geom_segment(aes(x = -0.7, y = -0.25, xend = -0.87, yend = -0.08), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(z[3222]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = -0.68, xmax = -0.68, 
                    ymin = -0.29, ymax = -0.29) +
  geom_segment(aes(x = -0.4, y = -0.25, xend = -0.23, yend = -0.08), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(z[2212]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = -0.38, xmax = -0.38, 
                    ymin = -0.29, ymax = -0.29) +
  geom_segment(aes(x = 1, y = -0.25, xend = 0.88, yend = -0.08), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(z[1212]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = 1, xmax = 1, 
                    ymin = -0.29, ymax = -0.29) +
  geom_segment(aes(x = 1.43, y = -0.25, xend = 1.55, yend = -0.08), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(z[3212]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = 1.43, xmax = 1.43, 
                    ymin = -0.29, ymax = -0.29) +
  geom_segment(aes(x = 2.12, y = -0.25, xend = 2, yend = -0.08), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(z[4222]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = 2.12, xmax = 2.12, 
                    ymin = -0.29, ymax = -0.29) +
  geom_segment(aes(x = 2.78, y = -0.25, xend = 2.9, yend = -0.08), arrow = arrow(length = unit(0.2, "cm"))) +
  annotation_custom(grob = textGrob(bquote(z[4212]),
                                    gp = gpar(fontsize = 20)), 
                    xmin = 2.78, xmax = 2.78, 
                    ymin = -0.29, ymax = -0.29)

linearPlotAFROC <- ggplotGrob(linearPlotAFROC)
linearPlotAFROC$layout$clip[linearPlotAFROC$layout$name=="panel"] <- "off"
grid.draw(linearPlotAFROC)
