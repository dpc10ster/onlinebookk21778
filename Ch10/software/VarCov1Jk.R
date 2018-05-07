VarCov1Jk <- function (JackFoMMatrix) 
{ 
  I <- dim(JackFoMMatrix)[1]
  K <- dim(JackFoMMatrix)[2]
  Cov <- array(dim = c(I, I))
  
  for (i in 1:I){
    for (ip in 1:I){
      Cov[i, ip] <- cov(JackFoMMatrix[i,], JackFoMMatrix[ip,])          
    }
  }  
  
  Var <- 0
  count <- 0
  for (i in 1:I){    
    Var <- Var + Cov[i, i] 
    count <- count + 1
  }
  Var <- Var / count 
  Var <- Var * (K-1)^2/K
  
  Cov1 <- 0
  count <- 0
  for (i in 1:I){    
    for (ip in 1:I){
      if (ip != i){
        Cov1 <- Cov1 + Cov[i, ip] 
        count <- count + 1
      }
    }
  }  
  Cov1 <- Cov1 / count 
  Cov1 <- Cov1 * (K-1)^2/K
  return (list (    
    Var = Var,
    Cov1 = Cov1
  ))  
}

