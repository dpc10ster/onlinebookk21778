rm( list = ls()) # mainRocPdfsWithCutoffs.R
require(ggplot2)
require(grid)

mu <- 1.5;sigma <- 1.5

z1 <- seq(-3, 4, by = 0.01)
z2 <- seq(-3, 6, by = 0.01)

Pdf1 <- dnorm(z1)
Pdf2 <- dnorm(z2, mu, sd = sigma)

df <- data.frame(z = c(z1, z2), pdfs = c(Pdf1, Pdf2), 
                 truth = c(rep('non-diseased', length(Pdf1)), rep('diseased', length(Pdf2))), 
                 stringsAsFactors = FALSE)

cut_point <- data.frame(z = c(-2.0, -0.5, 1, 2.5))

rocPdfs <- ggplot(df, aes(x = z, y = pdfs, color = truth)) + 
  geom_line(size = 2) + 
  scale_colour_manual(values=c("darkgrey","black")) + 
  theme(legend.title = element_blank(), legend.position = c(0.85, 0.95), 
        legend.text = element_text(size=25, face = "bold"), 
          axis.title.x = element_text(hjust = 0.8, size = 30,face="bold"),
        axis.title.y = element_text(size = 25,face="bold")) +
  geom_vline(data = cut_point, aes(xintercept = z), linetype = "dotted", size = 1.5) +
  annotation_custom(grob = textGrob(bquote(italic("O")),
                                    gp = gpar(fontsize = 32)), 
                    xmin = -3.2, xmax = -3.2, # adjust the position of "O"
                    ymin = -0.0, ymax = -0.01) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))


for (i in 1 : length(cut_point$z)){
  rocPdfs <- rocPdfs +
    annotation_custom(grob = textGrob(bquote(zeta[.(i)]),gp = gpar(fontsize = 20)),
                      xmin = cut_point$z[i], xmax = cut_point$z[i],
                      ymin = -0.01, ymax = -0.01)
}

gt <- ggplot_gtable(ggplot_build(rocPdfs))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

