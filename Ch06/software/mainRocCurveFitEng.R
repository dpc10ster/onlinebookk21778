#mainRocCurveFitEng.R
rm(list = ls())
source("rocY.R")
library(ggplot2);require("grid")
# the famous a and b parameters of the 
# Dorfman and Alf binormal model, 
# as yielded by the Eng JAVA code
a <- 1.3204; b <- 0.6075

FPF <- seq(0.0, 1, 0.01)
curveData <- data.frame(
  FPF = FPF, 
  TPF = rocY(FPF, a, b))
# observed operating points
# Table 6.1
FPF <- c(0.017, 0.050, 0.183, 0.5) 
TPF <- c(0.440, 0.680, 0.780, 0.900)
pointsData <- data.frame(
  FPF = FPF, TPF = TPF)
p <- ggplot(mapping = aes(x = FPF, y = TPF)) +
  theme(
    axis.title.y = element_text(size = 25,face="bold"),
    axis.title.x = element_text(size = 30,face="bold")) +
  geom_line(data = curveData, size = 2) +
  geom_point(data = pointsData, size = 5) +
  annotation_custom(
    grob = textGrob(bquote(italic("O")),
                    gp = gpar(fontsize = 32)), 
    xmin = -0.03, xmax = -0.03, # adjust the position of "O"
    ymin = -0.03, ymax = -0.03) +
  scale_x_continuous(
    expand = c(0, 0), 
    breaks = c(0.25, 0.5, 0.75, 1)) + 
  scale_y_continuous(
    expand = c(0, 0), 
    breaks = c(0.25, 0.5, 0.75, 1))
print(p)

p <- ggplotGrob(p)
p$layout$clip[p$layout$name=="panel"] <- "off"
grid.draw(p)

