Cov1Cov2Cov3 <- function (PseudovalueMatrix) 
{ 
  I <- dim(PseudovalueMatrix)[1]
  J <- dim(PseudovalueMatrix)[2]
  Covariance <- array(dim = c(I, I, J, J))
  
  for (i in 1:I){
    for (ip in 1:I){
      for ( j in 1:J) {   
        for ( jp in 1:J) {           
          Covariance[i, ip, j, jp] <- cov(PseudovalueMatrix[i, j, ], PseudovalueMatrix[ip, jp, ])          
        }
      }
    }
  }  
  
  Var <- 0
  count <- 0
  for (i in 1:I){    
    for (j in 1:J) {      
      Var <- Var + Covariance[i, i, j, j] 
      count <- count + 1
    }
  }
  Var <- Var / count 

  Vari <- array(0, dim = I)
  for (i in 1:I){    
    count <- 0
    for (j in 1:J) {      
      Vari[i] <- Vari[i] + Covariance[i, i, j, j] 
      count <- count + 1
    }
    Vari[i] <- Vari[i] / count   
  }

  
  Cov1 <- 0
  count <- 0
  for (i in 1:I){    
    for (ip in 1:I){
      for (j in 1:J) {      
        if (ip != i){
          Cov1 <- Cov1 + Covariance[i, ip, j, j] 
          count <- count + 1
        }
      }
    }
  }  
  Cov1 <- Cov1 / count 
  
  Cov2 <- 0
  count <- 0
  for (i in 1:I){    
    for (j in 1:J) {      
      for (jp in 1:J){
        if (j != jp){
          Cov2 <- Cov2 + Covariance[i, i, j, jp] 
          count <- count + 1
        }
      }
    }
  }  
  #Cov2 <- Cov2 / (I*J*(J-1)) # OK, DPC
  Cov2 <- Cov2 / count 
  
  # this is needed for single modality stats
  Cov2i <- array(0, dim = I)
   for (i in 1:I){    
     count <- 0
     for (j in 1:J) {      
      for (jp in 1:J){
        if (j != jp){
          Cov2i[i] <- Cov2i[i] + Covariance[i, i, j, jp] 
          count <- count + 1
        }
      }
    }
    Cov2i[i] <- Cov2i[i] / count 
   }  
  
  Cov3 <- 0
  count <- 0
  for (i in 1:I){
    for (ip in 1:I){
      if (i != ip){
        for (j in 1:J) {
          for (jp in 1:J){
            if (j != jp) {
              Cov3 <- Cov3 + Covariance[i, ip, j, jp] 
              count <- count + 1
            }
          }
        }
      }
    }
  }
  
  #Cov3 <- Cov3 / (I*(I-1)*J*(J-1)) # not OK; general advice; better to let computer do the thinking
  Cov3 <- Cov3 / count
  
  return (list (    
    Var = Var,
    Cov1 = Cov1,
    Cov2 = Cov2,
    Cov3 = Cov3,
    Vari = Vari,
    Cov2i = Cov2i
  ))  
}

