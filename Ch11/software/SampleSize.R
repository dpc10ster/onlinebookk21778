SampleSize <- function(J, KStar, VarTR, Var, Cov1, Cov2, Cov3, d,  
                                 alpha = 0.05, TargetPower = 0.8)
{
  K <- 20
  Power <- 0
  while(Power <= TargetPower){
    if (K > 2000){
      break
    }
    K <- K + 1
    Delta <- J*d^2/2/(VarTR+(KStar/K)*(Var-Cov1+(J-1)*max(Cov2-Cov3,0)))
    ddf <- (J-1)*(VarTR+(KStar/K)*(Var-Cov1+(J-1)*max(Cov2-Cov3,0)))^2/
      ((VarTR+(KStar/K)*(Var-Cov1-max(Cov2-Cov3,0)))^2)
    FCrit <- qf(1 - alpha, 1, ddf)
    Power <- 1- pf(FCrit, 1, ddf, ncp = Delta)
  }
  return(list(K = K, Power = Power))
}
