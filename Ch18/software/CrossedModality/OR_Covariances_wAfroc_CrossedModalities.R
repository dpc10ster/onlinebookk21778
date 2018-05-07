OR_Covariances_wAfroc_CrossedModalities <- function (JACK,fom,NL,LL,Nk,Wk,AvgIndx){
 
  if (!((fom == toupper("wAfroc")) || (fom == toupper("wAfroc1")))) stop("incorrect flag")
  
  I1 <- length(NL[,1,1,1,1])
  I2 <- length(NL[1,,1,1,1])
  J <- length(NL[1,1,,1,1])
  K2 <- length(LL[1,1,1,,1])
  K1 <- length(NL[1,1,1,,1]) - K2
  K <- K1 + K2
  
  
  if (JACK) {
    FOM_i1i2jk <- array(dim = c(I1,I2,J,K))
    for (k in 1:K1) {
      NL_jk <- NL[,,,-k,]
      LL_jk <- LL
      Nk_jk <- Nk
      Wk_jk <- Wk
      if (fom == toupper("wAfroc")) FOM_i1i2jk[,,,k] <- TrapezoidalAreawAfroc (NL_jk,LL_jk,Nk_jk,Wk_jk) 
      if (fom == toupper("wAfroc1")) FOM_i1i2jk[,,,k] <- TrapezoidalAreawAfroc1 (NL_jk,LL_jk,Nk_jk,Wk_jk) 
    }
    for (k in (K1+1):(K1+K2)) {
      NL_jk <- NL[,,,-k,]
      LL_jk <- LL[,,,-(k-K1),]
      Nk_jk <- Nk[-(k-K1)]
      Wk_jk <- Wk[-(k-K1),]
      if (fom == toupper("wAfroc")) FOM_i1i2jk[,,,k] <- TrapezoidalAreawAfroc (NL_jk,LL_jk,Nk_jk,Wk_jk) 
      if (fom == toupper("wAfroc1")) FOM_i1i2jk[,,,k] <- TrapezoidalAreawAfroc1 (NL_jk,LL_jk,Nk_jk,Wk_jk) 
    }    
    if (AvgIndx == 1) {
      FOM_i2jk <- array(dim = c(I2,J,K))
      for (i2 in 1:I2) {
        for (j in 1:J) {        
          for (k in 1:K) {        
            FOM_i2jk[i2,j,k] <- mean(FOM_i1i2jk[,i2,j,k])
          }
        }
      }
    } else {
      FOM_i1jk <- array(dim = c(I1,J,K))
      for (i1 in 1:I1) {
        for (j in 1:J) {        
          for (k in 1:K) {        
            FOM_i1jk[i1,j,k] <- mean(FOM_i1i2jk[i1,,j,k])
          }
        }
      }      
    }
    if (AvgIndx == 1) Cov <- Cov1Cov2Cov3(FOM_i2jk) else Cov <- Cov1Cov2Cov3(FOM_i1jk)
    Var = Cov$Var*(K-1)^2/K # see paper by Efron and Stein
    Cov1 = Cov$Cov1*(K-1)^2/K
    Cov2 = Cov$Cov2*(K-1)^2/K
    Cov3 = Cov$Cov3*(K-1)^2/K    
    Vari = Cov$Vari*(K-1)^2/K 
    Cov2i = Cov$Cov2i*(K-1)^2/K        
  } else {
    stop("Not implemented yet")
    FOM_b <- array(dim = c(I, J, B))    
    for (b in 1:B) {
      k_b <- ceiling( runif( K ) * K )
      ratings_b <- ratings[, , -k_b]
      truth_b <- truth[-k_b]
      if (fom == toupper("wAfroc")) FOM_b[, , b] <- TrapezoidalAreawAfroc(NL_b,LL_b,Nk_b,Wk_b) 
      if (fom == toupper("wAfroc1")) FOM_b[, , b] <- TrapezoidalAreawAfroc1(NL_b,LL_b,Nk_b,Wk_b) 
    }
    Cov <- Cov2Cov3(FOM_b)
    Var = Cov$Var
    Cov1 = Cov$Cov1
    Cov2 = Cov$Cov2
    Cov3 = Cov$Cov3
    Vari = Cov$Vari
    Cov2i = Cov$Cov2i    
  }
  
  if (fom == toupper("wAfroc")) FOM_i1i2j <- TrapezoidalAreawAfroc (NL,LL,Nk,Wk)
  if (fom == toupper("wAfroc1")) FOM_i1i2j <- TrapezoidalAreawAfroc1 (NL,LL,Nk,Wk)

  if (AvgIndx == 1) {
    FOM_i2j <- array(dim=c(I2,J))
    for (i2 in 1:I2) {
      for (j in 1:J) {                
        FOM_i2j[i2,j] <- mean(FOM_i1i2j[,i2,j])
      }
    }
  } else {
    FOM_i1j <- array(dim=c(I1,J))
    for (i1 in 1:I1) {
      for (j in 1:J) {                
        FOM_i1j[i1,j] <- mean(FOM_i1i2j[i1,,j])
      }
    }
  }
  
  if (AvgIndx == 1) {
    FOM <- FOM_i2j
  } else {
    FOM <- FOM_i1j 
  } 
    
  return (list(
    FOM = FOM,
    Var = Var,
    Cov1 = Cov1,
    Cov2 = Cov2,
    Cov3 = Cov3,
    Vari = Vari,
    Cov2i = Cov2i
  ))
}