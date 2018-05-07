# mainConfidenceIntervals.R
rm( list = ls() )
options(digits=6)
seed <- 100
set.seed(seed)
alpha <- 0.05
K1 <- 99
K2 <- 111
mu <- 5
zeta <- mu/2
cat("alpha = ", alpha, 
    "\nK1 = ", K1, 
    "\nK2 = ", K2, 
    "\nmu = ", mu, 
    "\nzeta = ", zeta, "\n")
z1 <- rnorm(K1)
z2 <- rnorm(K2) + mu
nTN <- length(z1[z1 < zeta])
nTP <- length(z2[z2 >= zeta])
Sp <- nTN/K1
Se <- nTP/K2
cat("Specificity = ", Sp, 
    "\nSensitivity = ", Se, "\n")

# Approx binomial tests
cat("approx 95% CI on Sp = ", 
    -abs(qnorm(alpha/2))*sqrt(Sp*(1-Sp)/K1)+Sp, 
    +abs(qnorm(alpha/2))*sqrt(Sp*(1-Sp)/K1)+Sp,"\n")

# Exact binomial test
ret <- binom.test(nTN, K1, p = nTN/K1)
cat("Exact 95% CI on Sp = ", 
    as.numeric(ret$conf.int),"\n")

# Approx binomial tests
cat("approx 95% CI on Se = ", 
    -abs(qnorm(alpha/2))*sqrt(Se*(1-Se)/K2)+Se, 
    +abs(qnorm(alpha/2))*sqrt(Se*(1-Se)/K2)+Se,"\n")

# Exact binomial test
ret <- binom.test(nTP, K2, p = nTP/K2)
cat("Exact 95% CI on Sp = ", 
    as.numeric(ret$conf.int),"\n")