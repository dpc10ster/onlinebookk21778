# mainPdfs.R
rm( list = ls())
require(ggplot2)
source("rocY.R")

aArray <- c(0.7, 0.7, 1.5, 2)
bArray <- c(0.5, 1.5, 0.5, 0.5)
z1 <- seq(-5, 3, by = 0.01)
z2 <- seq(-5, 7, by = 0.01)
FPF <- seq(0.0, 1, 0.01)

Pdf1 <- dnorm(z1)

for (i in 1:length(aArray))
{
  a <- aArray[i]
  b <- bArray[i]
  TPF <- rocY(FPF, a, b)
  TPF <- rocY(FPF, a, b)
  rocPlot <- data.frame(FPF = FPF, TPF = TPF)
  plotRoc <- ggplot(rocPlot, aes(x = FPF, y = TPF)) + 
    geom_line(size = 2)  + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"))  +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  print(plotRoc)
  
  Pdf2 <- dnorm(z2, a/b, sd = 1/b)
  df <- data.frame(z = c(z1, z2), pdfs = c(Pdf1, Pdf2), 
                   truth = c(rep('non-diseased', length(Pdf1)), 
                             rep('diseased', length(Pdf2))))
  
  rocPdfs <- ggplot(df, aes(x = z, y = pdfs
                            , color = truth)) + 
    geom_line(size = 2) + 
    scale_colour_manual(values=c("darkgrey","black")) + 
    theme(legend.title = element_blank(), legend.text=element_text(size=20,face="bold"), 
          legend.position = c(0.7, 0.9),
          legend.direction = "horizontal") + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"),
          legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(0.75, "inch"))  +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  print(rocPdfs)
  
  cat("a = ", a, ", b = ", b, "\n")
  break
}