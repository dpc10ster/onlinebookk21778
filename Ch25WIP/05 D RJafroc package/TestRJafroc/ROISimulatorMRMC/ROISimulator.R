
ROISimulator <- function( I, J, K, Q, mu, Deltamu, Var_R, Var_C, Var_mR, 
                              Var_mC, Var_RC, Var_e, rho_C, rho_RC, rho_mC, rho_e) 
{
  
  mu_it <- array(0, dim = c(I, 2))
  mu_it[,2] <- mu + Deltamu
  
  qk2 <- round(runif (K[2] , min = 1 , max = Q))
  
  isDiseased <- array(0 ,dim = c(K[2] , Q))
  for (k in 1:K[2]){
    isDiseased[k, sample(c(1:Q) , qk2[k] , replace = FALSE)] <- 1
  }
  
  Rjt <- rnorm( 2*J, sd = sqrt(Var_R) )
  dim(Rjt) <- c(J,2)
  
  Cktt <- rnorm( 2*max(K), sd = sqrt(Var_C * rho_C) )
  dim(Cktt) <- c(max(K),2)
  
  CLkttlss <- rnorm( 2*max(K)*Q, sd = sqrt(Var_C * (1 - rho_C)))
  dim(CLkttlss) <- c(max(K), 2, Q)
  
  mRijt <- rnorm( I*2*J, sd = sqrt(Var_mR) )
  dim(mRijt) <- c(I, J, 2)
  
  mCiktt <- rnorm( I*2*max(K), sd = sqrt(Var_mC * rho_mC) )
  dim(mCiktt) <- c(I, max(K), 2)
  
  mCLikttlss <- rnorm( I*2*max(K)*Q, sd = sqrt(Var_mC * (1 - rho_mC)) )
  dim(mCLikttlss) <- c(I, max(K), 2, Q)
  
  RCjktt <- rnorm( J*2*max(K), sd = sqrt(Var_RC * rho_RC) )
  dim(RCjktt) <- c(J, max(K), 2)
  
  RCLjkttlss <- rnorm( J*2*max(K)*Q, sd = sqrt(Var_RC * (1 - rho_RC)) )
  dim(RCLjkttlss) <- c(J, max(K), 2, Q)
  
  eijktt <- rnorm( I*J*2*max(K), sd = sqrt(Var_e * rho_e) )
  dim(eijktt) <- c(I, J, max(K), 2)
  
  eLijkttlss <- rnorm( I*J*2*max(K)*Q, sd = sqrt(Var_e * (1 - rho_e)) )
  dim(eLijkttlss) <- c(I, J, max(K), 2, Q)
  
  Rijkttlss <- array(dim=c(I, J,  max(K), 2, Q))
  for (i in 1:I) {
    for (j in 1:J) {
      for (t in 1:2) {
        if(t == 1){
          for (k in 1:K[t]) {
            for (r in 1:Q){ 
              Rijkttlss[i,j,k,t,r] <- (mu_it[i, t] + Rjt[j, t] + Cktt[k, t]  + CLkttlss[k, t, r]
                                         + mRijt[i, j, t] + mCiktt[i, k, t] + mCLikttlss[i, k, t, r] 
                                         + RCjktt[j, k, t] + RCLjkttlss[j, k, t, r] + eijktt[i,j,k,t] + eLijkttlss[i, j, k, t, r])
            }            
          }
        }else{
          for (k in 1:K[t]) {
            for (r in 1:Q){ 
              if (isDiseased[k, r] == 0) {
                Rijkttlss[i,j,k,t,r] <- (mu_it[i, 1] + Rjt[j, t] + Cktt[k, t]  + CLkttlss[k, t, r]
                                           + mRijt[i, j, t] + mCiktt[i, k, t] + mCLikttlss[i, k, t, r] 
                                           + RCjktt[j, k, t] + RCLjkttlss[j, k, t, r] + eijktt[i,j,k,t] + eLijkttlss[i, j, k, t, r])
              } else {
                Rijkttlss[i,j,k,t,r] <- (mu_it[i, t] + Rjt[j, t] + Cktt[k, t]  + CLkttlss[k, t, r]
                                           + mRijt[i, j, t] + mCiktt[i, k, t] + mCLikttlss[i, k, t, r] 
                                           + RCjktt[j, k, t] + RCLjkttlss[j, k, t, r] + eijktt[i,j,k,t] + eLijkttlss[i, j, k, t, r])
              }
            }            
          }
        }
      }  
    }
  }
  
  return( list (Rijkttlss = Rijkttlss,
                isDiseased = isDiseased
  ))
}
