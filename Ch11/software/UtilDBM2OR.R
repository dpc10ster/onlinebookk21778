UtilDBM2OR <- function(I, J, K, MSC,MSTC,MSRC,MSTRC) {
  
  Var <- (MSC+(I-1)*MSTC+(J-1)*MSRC+(I-1)*(J-1)*MSTRC)/I/J/K
  Cov1 <- (MSC-MSTC+(J-1)*(MSRC-MSTRC))/I/J/K
  Var <- (MSC-MSTC+(I-1)*(MSTC-MSTRC))/I/J/K
  Var <- (MSC-MSTC-MSRC+MSTRC)/I/J/K
  return (list (    
    Var = Var,
    Cov1 = Cov1,
    Cov2 = Cov2,
    Cov3 = Cov3    
  ))  
}
  