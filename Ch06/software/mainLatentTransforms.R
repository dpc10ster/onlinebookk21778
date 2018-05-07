# mainLatentTransforms.R
# shows that monontone transformations have no effect on 
# AUC even though the pdfs look non-gaussian
# common statistician misconception about ROC analysis

rm( list = ls( all = TRUE ) )
library(RJafroc)
library(ggplot2)
source( 'TrapezoidalArea.R' )
options(digits = 7)

Y <- function(z,mu1,mu2,sigma1,sigma2,f) {
  y <- (1-f)*pnorm((z-mu1)/sigma1)*100
  +f*pnorm((z-mu2)/sigma2)*100
  return( y )
}

fArray <- c(0.1,0.5,0.9)
seedArray <- c(10,11,12)
for (i in 1:3) { 
  f <- fArray[i];seed <- seedArray[i];set.seed(seed)
  K1 <- 900;K2 <- 1000
  mu1 <- 30;sigma1 <- 7;mu2 <- 55;sigma2 <- 7
  z1 <- rnorm(K1, mean = mu1, sd = sigma1)
  z1[z1>100] <- 100;z1[z1<0] <- 0
  z2 <- rnorm(K2, mean = mu2, sd = sigma2)
  z2[z2>100] <- 100;z2[z2<0] <- 0
  AUC1 <- TrapezoidalArea(z1, z2)
  Gaussians <- c(z1, z2)
  hist1 <- data.frame(x=Gaussians)
  hist.1 <-  ggplot(data = hist1, mapping = aes(x = x)) +
    geom_histogram(
      binwidth = 1, 
      color = "black", 
      fill="grey") + 
    xlab(label = "Original Rating") + 
    ggtitle(label = "Gaussians")
  print(hist.1)
  
  z <- seq(0.0, 100, 0.1)
  curveData <- data.frame(
    x = z, 
    z =  Y(z,mu1,mu2,sigma1,sigma2,f))
  plot3 <- ggplot(mapping = aes(x = x, y = z)) + 
    geom_line(data = curveData) +
    xlab(label = "Original Rating") +
    ylab(label = "Transformed Rating") + 
    ggtitle(label = "Monotone Transformation")
  print(plot3)
  
  y <- Y(c(z1, z2),mu1,mu2,sigma1,sigma2,f)
  y1 <- y[1:K1];y2 <- y[(K1+1):(K1+K2)]
  AUC2 <- TrapezoidalArea( y1, y2)
  hist2 <- data.frame(x=y)
  hist.2 <-  ggplot(data = hist2, mapping = aes(x = x)) +
    geom_histogram(
      binwidth = 1, 
      color = "black", 
      fill="grey") + 
    xlab(label = "Transformed Rating") + 
    ggtitle(label = "Latent Gaussians")
  print(hist.2)
  cat("seed =", seed, 
      "\nf =", f, 
      "\nAUC of actual Gaussians =", AUC1, 
      "\nAUC of latent Gaussians =", AUC2, "\n")
}