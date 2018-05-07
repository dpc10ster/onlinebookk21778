VarCov1Bs <- function (FomBs) 
{ 
  I <- dim(FomBs)[1]
  Covariance <- array(dim = c(I, I))
  
  for (i in 1:I){
    for (ip in 1:I){
      Covariance[i, ip] <- cov(FomBs[i, ], FomBs[ip, ])          
    }
  }  
  #CovarianceFirstPrinciples(FomBs)
  Var <- 0
  count <- 0
  for (i in 1:I){    
    Var <- Var + Covariance[i, i] 
    count <- count + 1
  }
  Var <- Var / count 
  
  Cov1 <- 0
  count <- 0
  for (i in 1:I){    
    for (ip in 1:I){
      if (ip != i){
        Cov1 <- Cov1 + Covariance[i, ip] 
        count <- count + 1
      }
    }
  }  
  Cov1 <- Cov1 / count 
  
  return (list (    
    Var = Var,
    Cov1 = Cov1
  ))  
}

