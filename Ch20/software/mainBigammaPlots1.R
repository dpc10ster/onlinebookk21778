rm(list = ls()) #mainBigammaPlots.R
require(stats)
require(ggplot2)

aucIntegrand <- function (FPF, r, lambda) 
{
  y <- 1 - pgamma(qgamma(1-FPF, r), r, scale = 1/lambda)
  return(y)  
}


biGammaRocY <- function (x, r, lambda) {
  y <- 1 - pgamma(x, r, scale = 1/lambda)
  return(y)
}

rArray <- c(1,4.391,5,10);lambdaArray <- c(1,0.439,0.3,0.1)
#r <-  1; lambda <-  1 # I made it up
#r <-  4.391; lambda <-  0.439 # from Dorfman paper
#r <-  5; lambda <-  0.3 # I made it up
#r <-  10; lambda <-  0.1 # I made it up
# left limit below gets the upper end of the ROC curve; the right limit gets the lower corner
for (i in 1:length(rArray))
{
  r <- rArray[i];lambda <- lambdaArray[i]
  x <- seq(0, r/lambda + 10*sqrt(r/lambda^2), by = 0.01)
  FPF <- biGammaRocY(x, r, 1)
  TPF <- biGammaRocY(x, r, lambda)
  
  rocPlot <- data.frame(FPF = FPF, TPF = TPF)
  plotRoc <- ggplot(rocPlot, aes(x = FPF, y = TPF)) + geom_line() + 
    geom_line(data = rocPlot, size = 4) + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"))  +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  print(plotRoc)
  
  x <- seq(0.01, r/lambda + 3*sqrt(r/lambda^2), by = 0.01) # taking mean + 3 times the standard deviation
  
  Pdf1 <- dgamma(x, r)
  Pdf2 <- dgamma(x, r, scale = 1/lambda)
  
  df <- data.frame(x = c(x, x), pdf = c(Pdf1, Pdf2), 
                   truth = c(rep('non-diseased', length(Pdf1)), 
                             rep('diseased', length(Pdf2))))
  
  bigammapdfs <- ggplot(df, aes(x = x, y = pdf, color = truth)) + 
    geom_line(size = 2) +
    scale_color_manual(values = c("black", "darkgrey")) + 
    theme(legend.title = element_blank(), legend.position = c(0.75, 0.95),
          legend.key.size = unit(0.5, "inch"), legend.key.height = unit(0.5, "inch"),  legend.text=element_text(size=20, face = "bold"),
          legend.direction = "horizontal") +
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  
  print(bigammapdfs)
  
  AUC <- integrate(aucIntegrand,0,1, r = r, lambda = lambda)
  cat("r = ", r, ", lambda = ", lambda, ", AUC = ", AUC$value, "\n")
  
}