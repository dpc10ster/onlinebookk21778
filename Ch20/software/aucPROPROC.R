aucPROPROC <- function (c1, da){
  # Metz and Pan Journal of Mathematical Psychology 43, 1?33 (1999)
  rho2 <- -(1-c1^2)/(1+c1^2)
  corr <- diag(2)
  corr[lower.tri(corr)] <- rho2
  corr[upper.tri(corr)] <- rho2
  lower <- rep(-Inf,2)
  upper <- c(-da/sqrt(2),0)
  mean <- rep(0,2)
  aucProproc <- pnorm(da/sqrt(2))+2*pmvnorm(lower, upper, mean, corr)  #Eqn. 36 Metz and Pan
  return (aucProproc)
}