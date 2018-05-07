rm(list = ls()) #MainDBMHBrief2.R  # freeze lines ##Fig. 9.5(a,b))
library(RJafroc);require("ggplot2");require("grid")
fileName <- "CXRinvisible3-20mm.xlsx"
frocData <- DfReadDataFile(fileName, format = "JAFROC", renumber = "TRUE")
rocData <- DfFroc2Roc(frocData)

plotM <- list(1, 2)
plotR <- list(c(1:5), c(1:5))
plot12Avg <- PlotEmpiricalOperatingCharacteristics(rocData, trts = plotM, rdrs = plotR)
p <- plot12Avg$Plot + scale_color_manual(values = c("black","darkgrey"))  +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(1,0), legend.direction = "horizontal",
        legend.text = element_text(size = 20, face = "bold"),legend.key.size = unit(2.5, "lines")) +
  annotation_custom(grob = textGrob(bquote(italic("O")),gp = gpar(fontsize = 32)), 
                    xmin = -0.03, xmax = -0.03, # adjust the position of "O"
                    ymin = -0.03, ymax = -0.03) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0,1), x = c(0,1))  
p$layers[[1]]$aes_params$size <- 2 # line
p$layers[[2]]$aes_params$size <- 5 # points

p <- ggplotGrob(p)
p$layout$clip[p$layout$name=="panel"] <- "off"
grid.draw(p)

plotM <- list(1, 3)
plotR <- list(c(1:5), c(1:5))
plot13Avg <- PlotEmpiricalOperatingCharacteristics(rocData, trts = plotM, rdrs = plotR)
p <- plot13Avg$Plot  + scale_color_manual(values = c("black","darkgrey"))  +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(1,0), legend.direction = "horizontal",
        legend.text = element_text(size = 20, face = "bold"),legend.key.size = unit(2.5, "lines")) +
annotation_custom(grob = textGrob(bquote(italic("O")),gp = gpar(fontsize = 32)), 
                  xmin = -0.03, xmax = -0.03, # adjust the position of "O"
                  ymin = -0.03, ymax = -0.03) +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) + 
  scale_y_continuous(expand = c(0, 0), breaks = c(0.25, 0.5, 0.75, 1)) +
  coord_cartesian(ylim = c(0,1), x = c(0,1))  
p$layers[[1]]$aes_params$size <- 2 # line
p$layers[[2]]$aes_params$size <- 5 # points

p <- ggplotGrob(p)
p$layout$clip[p$layout$name=="panel"] <- "off"
grid.draw(p)

