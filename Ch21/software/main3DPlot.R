rm(list = ls())
library(plotly)
library(plot3D)
library(mvtnorm)

muX <- 1.5
muY <- 2.3
sigmaX <- 1.5
sigmaY <- 2.6
rho1 <- 0.3
rho2 <- 0.9
SIGMA1 <- rbind(c(1, rho1), c(rho1, 1))
SIGMA2 <- rbind(c(sigmaX^2, sigmaX * sigmaY * rho2), c(sigmaX * sigmaY * rho2, sigmaY^2))

upperX <- ceiling(muX + 3 * sigmaX)
upperY <- ceiling(muY + 3 * sigmaY)
x <- seq(-4, upperX, by = 0.1)
y <- seq(-4, upperX, by = 0.1)
z <- array(dim = c(length(x), length(y)))
for (ix in 1:length(x)){
  for (iy in 1:length(y)){
    z[ix, iy] <- max(dmvnorm(c(x[ix], y[iy]), sigma = SIGMA1), dmvnorm(c(x[ix], y[iy]), mean = c(muX, muY), sigma = SIGMA2))
  }
}

p <- plot_ly(x = x, y = y, z = z, type = "surface")
