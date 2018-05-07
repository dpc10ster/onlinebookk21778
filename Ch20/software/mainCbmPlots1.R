rm(list = ls()) #mainCbmPlots.R
library(ggplot2)
source("cbmFunctions.R")

FPF <- seq(0.0, 1, 0.001)
alphaArr <- c(0.2, 0.8);muArr <- c(1,3)

for (i in 1:2)
  for (j in 1:2) 
  {
    {
      alpha <- alphaArr[i]
      mu <- muArr[j]
      TPF <- CbmRocY(FPF, mu, alpha)
      
      rocPlot <- data.frame(FPF = FPF, TPF = TPF)
      plotRoc <- ggplot(rocPlot, aes(x = FPF, y = TPF)) + geom_line() + 
        geom_line(data = rocPlot, size = 2) + 
        theme(axis.title.y = element_text(size = 25,face="bold"),
              axis.title.x = element_text(size = 30,face="bold"))  +
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) 
      #print(plotRoc)
    
      if (i == 1) {
        z1 <- seq(-3, 3, by = 0.01)
        z2 <- seq(-3, mu + 3, by = 0.01)
      } else {
        z1 <- seq(-3, 3, by = 0.01)
        z2 <- seq(-3, mu + 3, by = 0.01)
      }
      
      Pdf1 <- dnorm(z1)
      Pdf2 <- (1 - alpha) * dnorm(z2) + alpha * dnorm(z2, mu)
      
      df <- data.frame(
        z = c(z1, z2), pdf = c(Pdf1, Pdf2), 
        truth = c(rep('non-diseased', length(Pdf1)), 
                  rep('diseased', length(Pdf2)))
      )
      
      cbmPdfs <- ggplot(df, aes(x = z, y = pdf, color = truth)) +
        geom_line(data = df, size = 2) +
        scale_colour_manual(values=c("black","darkgrey")) +
        theme(legend.title = element_blank(), legend.position = c(0.5, 0.1),
              legend.key.size = unit(0.5, "inch"), legend.key.height = unit(0.5, "inch"),  legend.text=element_text(size=20, face = "bold"),
              legend.direction = "horizontal") +
        theme(axis.title.y = element_text(size = 25,face="bold"),
              axis.title.x = element_text(size = 30,face="bold"))  +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))

      print(cbmPdfs)
      cat("mu = ", mu, ", alpha = ", alpha, "\n")
      next
    }
  }