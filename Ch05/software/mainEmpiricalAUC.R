rm( list = ls())#mainEmpiricalAUC.R
require("ggplot2");require("grid")
source("RocOperatingPointsFromRatingsTable.R")

RocDataTable = array(dim = c(2,4))
RocDataTable[1,]  <- c(30,19,8,3)
RocDataTable[2,]  <- c(5,11,12,22)

ret  <- RocOperatingPointsFromRatingsTable( RocDataTable[1,], RocDataTable[2,] )
FPF <- ret$FPF;TPF <- ret$TPF

ROC_Points <- data.frame(FPF = FPF, TPF = TPF)
ROC_Points <- rbind(c(0, 0), ROC_Points, c(1, 1)) # add the trivial points


shade <- data.frame(FPF = c(ROC_Points$FPF, 1), TPF = c(ROC_Points$TPF, 0))

p <- ggplot(ROC_Points, aes(x = FPF, y = TPF) ) + 
  geom_polygon(data = shade, fill = 'grey') + 
  geom_line(size = 2) + 
  geom_point(size = 5) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = expression(FPF)) + labs(y = expression(TPF)) + 
  coord_cartesian(ylim = c(0,1), x = c(0,1)) + xlab("FPF")+ ylab("TPF" ) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0,1), x = c(0,1)) + 
  annotation_custom(grob = textGrob(bquote(italic("O")),
                    gp = gpar(fontsize = 32)), 
                    xmin = -0.03, xmax = -0.03, # adjust the position of "O"
                    ymin = -0.03, ymax = -0.03) +
  theme(axis.text=element_text(size = 25,face="bold"), axis.title=element_text(size=30,face="bold"))

p <- ggplotGrob(p)
p$layout$clip[p$layout$name=="panel"] <- "off"
grid.draw(p)
