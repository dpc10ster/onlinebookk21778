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
        geom_line(data = rocPlot)
      print(plotRoc)
    
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
        geom_line(data = df)
      print(cbmPdfs)
      cat("mu = ", mu, ", alpha = ", alpha, "\n")
      next
    }
  }