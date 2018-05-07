# mainBivariateIntegrals.R
rm(list = ls()) 
library(mvtnorm)

cat("Integrals of the bivariate normal distribution\n")
mu  <- c(0,0)
rho <- 0.5;var1 <- 1.5;var2 <- 1.5;
COV <- matrix(
  c(var1, rho*sqrt(var1*var2), 
    rho*sqrt(var1*var2), var2),2)

ret <- pmvnorm(
  c(-Inf, -Inf), 
  c(Inf, Inf),
  mean=mu, 
  sigma = COV)
cat("Over the entire space = ", ret, "\n")
ret <- pmvnorm(
  c(-Inf, -Inf), 
  c(Inf, 0),
  mean=mu, 
  sigma = COV)
cat("Over the full space in one dimension\n") 
cat("and the -ve half space in other dimension = ", 
    ret, "\n")
ret <- pmvnorm(
  c(.3, .4), 
  c(.4, .5),
  mean=mu, 
  sigma = COV)
cat("Between specified ctff. values = ", 
    ret, "\n")
