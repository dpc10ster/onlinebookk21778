rm(list = ls()) # mainEmpRocPlot.R
library(ggplot2)
library(grid)
options(digits = 4)
K1 <- 60;K2 <- 50
FPF <- c(0, cumsum(rev(c(30, 19, 8, 2, 1))) / K1)
TPF <- c(0, cumsum(rev(c(5, 6, 5, 12, 22))) / K2)

ROCOp <- data.frame(FPF = FPF, TPF = TPF)
ROCPlot <- ggplot(data = ROCOp, mapping = aes(x = FPF, y = TPF)) + geom_line(size = 2) + 
  geom_point(size = 5) +  xlab("FPF")+ ylab("TPF" ) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black"), 
        axis.text = element_text(size = 25,face="bold"), 
        axis.title = element_text(size = 30,face="bold")) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0,1), x = c(0,1)) + 
  annotation_custom(grob = textGrob(bquote(italic("O")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = -0.03, xmax = -0.03, # adjust the position of "O"
                    ymin = -0.03, ymax = -0.03) +
  annotation_custom(grob = textGrob(bquote(italic(O[4])),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.06, xmax = 0.06, 
                    ymin = 0.43, ymax = 0.43) +
  annotation_custom(grob = textGrob(bquote(italic(O[3])),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.09, xmax = 0.09, 
                    ymin = 0.67, ymax = 0.67) +
  annotation_custom(grob = textGrob(bquote(italic(O[2])),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.18, xmax = 0.18, 
                    ymin = 0.82, ymax = 0.82) +
  annotation_custom(grob = textGrob(bquote(italic(O[1])),
                                    gp = gpar(fontsize = 32)), 
                    xmin = 0.49, xmax = 0.49, 
                    ymin = 0.93, ymax = 0.93)  

p <- ggplotGrob(ROCPlot)
p$layout$clip[p$layout$name=="panel"] <- "off"
grid.draw(p)