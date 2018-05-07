rm(list = ls()) #mainBigammaPlots.R
require(stats)
require(ggplot2)

aucIntegrand <- function (FPF, r, lambda) 
{
  y <- 1 - pgamma(qgamma(1-FPF, r), r, scale = 1/lambda)
  return(y)  
}


biGammaRocY <- function (xi, r, lambda) {
  y <- 1 - pgamma(xi, r, scale = 1/lambda)
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
  xi <- seq(0, r/lambda + 10*sqrt(r/lambda^2), by = 0.01)
  FPF <- biGammaRocY(xi, r, 1)
  TPF <- biGammaRocY(xi, r, lambda)
  
  rocPlot <- data.frame(FPF = FPF, TPF = TPF)
  plotRoc <- ggplot(rocPlot, aes(x = FPF, y = TPF)) + 
    geom_line() + 
    geom_line(data = rocPlot)
  print(plotRoc)
  
  xi <- seq(0.01, r/lambda + 3*sqrt(r/lambda^2), by = 0.01) # taking mean + 3 times the standard deviation
  
  Pdf1 <- dgamma(xi, r)
  Pdf2 <- dgamma(xi, r, scale = 1/lambda)
  
  df <- data.frame(xi = c(xi, xi), pdf = c(Pdf1, Pdf2), 
                   truth = c(rep('non-diseased', length(Pdf1)), 
                             rep('diseased', length(Pdf2))))
  
  bigammapdfs <- ggplot(df, aes(x = xi, y = pdf, color = truth)) + 
    geom_line() 
  print(bigammapdfs)
  
  AUC <- integrate(aucIntegrand,0,1, r = r, lambda = lambda)
  cat("r = ", r, ", lambda = ", lambda, ", AUC = ", AUC$value, "\n")
  
}