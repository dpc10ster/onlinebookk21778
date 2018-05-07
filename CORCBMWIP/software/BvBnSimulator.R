BvBnSimulator <- function(ax, bx, ay, by, rho1, rho2, K1, K2){
  sigma1 <- rbind(c(bx^2, rho1 * bx * by), c(rho1 * bx * by, by^2))
  sigma2 <- rbind(c(1, rho1), c(rho1, 1))
  z1 <- rmvnorm(K1, mean = c(0, 0), sigma = sigma1)
  z2 <- rmvnorm(K2, mean = c(ax, ay), sigma = sigma1)
  return(list(z1 = z1, 
              z2 = z2))
}