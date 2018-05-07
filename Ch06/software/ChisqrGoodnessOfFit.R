# parameters are ordered as a,b,zeta1,..., zetaLast
ChisqrGoodnessOfFit <- function (parameters,K1,K2)
{
  R <- length(K1);L <- length(parameters)
  
  zeta <- c(-Inf, parameters[3:L],Inf)
  a <- parameters[1]
  b <- parameters[2]
  k1 <- sum(K1);k2 <- sum(K2) # total number of non-diseased and diseased cases
  
  K1Exp <- rep(NA, length(K1))
  K2Exp <- rep(NA, length(K2))
  for (r in 1:R) {
    K1Exp[r] <- pnorm(zeta[r+1]) - pnorm(zeta[r])
    K2Exp[r] <- pnorm(b*zeta[r+1]-a) - pnorm(b*zeta[r]-a)
  }
  return(UtilGoodnessOfFit(rbind(K1, K2), rbind(K1Exp, K2Exp), 2))
}