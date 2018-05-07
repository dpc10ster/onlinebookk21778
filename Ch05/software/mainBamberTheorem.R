rm( list = ls()) #MainBamberTheorem.R
require("ggplot2");require("grid")
source("RocOperatingPoints.R")

RocCountsTable = array(dim = c(2,5))
RocCountsTable[1,]  <- c(30,19,8,2,1)
RocCountsTable[2,]  <- c(5,6,5,12,22)

ret  <- RocOperatingPoints( RocCountsTable[1,], 
                            RocCountsTable[2,] )
FPF <- ret$FPF;TPF <- ret$TPF

ROC_Points <- data.frame(FPF = FPF, TPF = TPF)
# Operating points without trivial ones
Op_Points <- ROC_Points
# add the trivial points
ROC_Points <- rbind(c(0, 0), ROC_Points, c(1, 1)) 

zeta1_x <- ROC_Points$FPF[4]
zeta1_y <- ROC_Points$TPF[4]
zeta1 <- data.frame(
  x = c(zeta1_x, zeta1_x, 1), 
  y = c(0, zeta1_y, zeta1_y))
i1 <- data.frame(x = zeta1_x, y = zeta1_y, label = "i + 1")
zeta2_x <- ROC_Points$FPF[5]
zeta2_y <- ROC_Points$TPF[5]
zeta2 <- data.frame(
  x = c(zeta2_x, zeta2_x, 1), 
  y = c(0, zeta2_y, zeta2_y))
i2 <- data.frame(x = zeta2_x, y = zeta2_y, label = "i")

shade <- data.frame(FPF = c(zeta2_x, zeta1_x, 1, 1), 
                    TPF = c(zeta2_y, zeta1_y, zeta1_y, zeta2_y))

p <- ggplot(ROC_Points, aes(x = FPF, y = TPF) ) + 
  geom_polygon(data = shade, fill = 'grey') + geom_line(size = 2) + 
  geom_point(data = Op_Points, size = 5) + 
  geom_line(data = zeta1, aes(x = x, y = y), linetype = 2, size = 2) + 
  geom_line(data = zeta2, aes(x = x, y = y), linetype = 2, size = 2) +  
  geom_text(data = i1, aes(x = x, y = y, label = label), size = 10, 
            fontface="bold", vjust = -0.5, hjust = 1) + 
  geom_text(data = i2, aes(x = x, y = y, label = label), size = 10, 
            fontface="bold", vjust = -0.5, hjust = 3) + 
  annotation_custom(grob = textGrob(bquote(bold(P(Z[1] >= zeta[i+1]))),
                                    gp = gpar(fontsize = 20)), 
                    xmin = zeta1_x - 0.08, xmax = zeta1_x - 0.08, 
                    ymin = 0.12, ymax = 0.12) +
  geom_segment(aes(x = zeta1_x - 0.05, y = 0.1, xend = zeta1_x, yend = 0), 
               arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  annotation_custom(grob = textGrob(bquote(bold(P(Z[1] >= zeta[i]))),
                                    gp = gpar(fontsize = 20)), 
                    xmin = zeta2_x + 0.08, xmax = zeta2_x + 0.08, 
                    ymin = 0.12, ymax = 0.12) +
  geom_segment(aes(x = zeta2_x + 0.05, y = 0.1, xend = zeta2_x, yend = 0), 
               arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  annotation_custom(grob = textGrob(bquote(bold(P(Z[2] >= zeta[i+1]))),
                                    gp = gpar(fontsize = 20)), 
                    xmin = 0.85, xmax = 0.95, 
                    ymin = zeta1_y - 0.12, ymax = zeta1_y - 0.12) +
  geom_segment(aes(x = 0.95, y = zeta1_y - 0.1, xend = 1, yend = zeta1_y), 
               arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  annotation_custom(grob = textGrob(bquote(bold(P(Z[2] >= zeta[i]))),
                                    gp = gpar(fontsize = 20)), 
                    xmin = 0.8, xmax = 0.8, 
                    ymin = zeta2_y + 0.05, ymax = zeta2_y + 0.02) +
  geom_segment(aes(x = 0.88, y = zeta2_y + 0.03, xend = 1, yend = zeta2_y), 
               arrow = arrow(angle = 20, length = unit(0.5, "cm"))) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotation_custom(grob = textGrob(bquote(italic("O")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = -0.03, xmax = -0.03, # adjust the position of "O"
                    ymin = -0.03, ymax = -0.03) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0,1), x = c(0,1)) +  
  xlab("FPF")+ ylab("TPF" ) + 
  theme(axis.text=element_text(size = 25,face="bold"), 
        axis.title=element_text(size=30,face="bold"))

# the following three lines are refreshing the plot, since "O" is out of the plotting area
p <- ggplotGrob(p)
p$layout$clip[p$layout$name=="panel"] <- "off"
grid.draw(p)
