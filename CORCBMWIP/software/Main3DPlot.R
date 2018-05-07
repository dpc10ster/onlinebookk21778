rm(list = ls())

library(plotly)
library(plot3D)
library(mvtnorm)

muX <- 3
muY <- 3
mu1 <- c(0, 0)
mu2 <- c(0, muY)
mu3 <- c(muX, 0)
mu4 <- c(muX, muY)

rho1 <- 0.3
rho4 <- 0.8
rho2 <- (rho1 + rho4) / 2
rho3 <- rho2
alphaX <- 0.5
alphaY <- 0.6

sigma1 <- rbind(c(1, rho1), c(rho1, 1))
sigma2 <- rbind(c(1, rho2), c(rho2, 1))
sigma3 <- rbind(c(1, rho3), c(rho3, 1))
sigma4 <- rbind(c(1, rho4), c(rho4, 1))

x <- seq(-3, 6, by = 0.05)
y <- x

M <- mesh(x, y)
X <- M$x
Y <- M$y
PDF <- (1 - alphaX) * (1 - alphaY) * dmvnorm(cbind(as.vector(X), as.vector(Y)), sigma = sigma1) + 
  (1 - alphaX) * alphaY * dmvnorm(cbind(as.vector(X), as.vector(Y)), mean = mu2, sigma = sigma2) + 
  alphaX * (1 - alphaY) * dmvnorm(cbind(as.vector(X), as.vector(Y)), mean = mu3, sigma = sigma3) + 
  alphaX * alphaY * dmvnorm(cbind(as.vector(X), as.vector(Y)), mean = mu4, sigma = sigma4)
dim(PDF) <- c(length(x), length(y))

scene <- list(camera = list(eye = list(x = -1.25, y = -2, z = 1.25)))

#p <- plot_ly(x = x, y = y, z = PDF, type = "surface") %>% layout(title = "PDF of Abnormal Cases", scene = scene)  %>% hide_colorbar()
p <- plot_ly(x = x, y = y, z = PDF, type = "surface") %>% hide_colorbar()
print(p)
