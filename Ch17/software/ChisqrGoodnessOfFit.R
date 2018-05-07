# parameters are ordered as a,b,zeta1,..., zetaLast
ChisqrGoodnessOfFit <- function (parameters,K1,K2)
{
  R <- length(K1);L <- length(parameters)
  
  zeta <- c(-Inf, parameters[3:L],Inf)
  a <- parameters[1]
  b <- parameters[2]
  k1 <- sum(K1);k2 <- sum(K2) # total number of non-diseased and diseased cases
  
  C2 <- 0
  for (r in 1:R) {
    K1Exp <- k1*(pnorm(zeta[r+1]) - pnorm(zeta[r]));if (K1[r] < 5){C2 <- NA;break}
    K2Exp <- k2*(pnorm(b*zeta[r+1]-a) - pnorm(b*zeta[r]-a));if (K2[r] < 5){C2 <- NA;break}    
    C2 <- C2 + (K1[r]-K1Exp)^2/K1Exp + (K2[r]-K2Exp)^2/K2Exp
  }
  
  if(!is.na(C2)) {
    df <- 2*R -(2+R-1) - 2 
    pVal <- 1- pchisq(C2,df)    
    return( list(
      C2 = C2,
      df = df,
      pVal = pVal
    ) )
  } else return (NA)  
}
