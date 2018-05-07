# MainBinaryRatings.R
rm( list = ls() ) # clear all variables
mu <- 1.5
zeta <- mu/2
seed <- 100 # initial value 100, then 101
set.seed(seed)
K1 <- 9
K2 <- 11
z1 <- rnorm(K1)
z2 <- rnorm(K2) + mu
nTN <- length(z1[z1 < zeta])
nTP <- length(z2[z2 >= zeta])
Sp <- nTN/K1
Se <- nTP/K2 
mu <- qnorm(Sp)+qnorm(Se)
cat("seed = ", seed, 
    "\nK1 = ", K1, 
    "\nK2 = ", K2, 
    "\nSpecificity = ", Sp, 
    "\nSensitivity = ", Se, 
    "\nEstimated mu = ", mu, "\n")

# # example 1
# pnorm(0.75)
# # example 2
# 2*qnorm(pnorm(zeta))
