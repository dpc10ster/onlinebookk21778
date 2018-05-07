TrapezoidalAreawAfroc1 <- function (NL, LL, Nk, Wk)
{
  MultipleModality <- FALSE
  NestedModality <- FALSE
  data <- dim(NL)
  if (length(data) == 5) {
    NestedModality <- TRUE
    I1 <- data[1]
    I2 <- data[2]
    J <- data[3]
    K1 <- length(NL[1,1,1,,1])
  } else if (length(data) == 4) {
    MultipleModality <- TRUE
    I <- data[1]
    J <- data[2]
    K1 <- length(NL[1,1,,1])
  } else if (length(data) == 3) {
    J <- data[1]
    K1 <- length(NL[1,,1])
  } else stop("Error in TrapezoidalAreaROC function")
  
  K2 <- length(Nk)
  K1 <- K1 - K2
  K <- K1 + K2
  
  if (NestedModality && !MultipleModality){
    theta <- array(0, dim = c(I1,I2,J))
    for (i1 in 1 : I1){
      for (i2 in 1 : I2){
        for (j in 1 : J) {
          for (k in 1:K){ #    
            X <- max(NL[i1,i2,j,k,]) # all cases are being used to get max NL rating
            theta[i1,i2,j] <- theta[i1,i2,j] + sum( (X < LL[i1,i2,j,,]) * Wk ) + 0.5 * sum( (X == LL[i1,i2,j,,]) * Wk)            
          }  
          theta[i1,i2,j] <- theta[i1,i2,j] / K / K2#
        }
      }
    }
  }
  
  if (!NestedModality && MultipleModality){
    theta <- array(0, dim = c(I,J))
    for (i in 1 : I){
      for (j in 1 : J) {
        for (k in 1:K) {#    
          X <- max(NL[i,j,k,]) # all cases are being used to get max NL rating
          theta[i,j] <- theta[i,j] + sum( (X < LL[i,j,,]) * Wk ) + 0.5 * sum( (X == LL[i,j,,]) * Wk)            
        }  
        theta[i,j] <- theta[i,j] / K / K2#
      }
    }
  }
  
  if (!NestedModality && !MultipleModality){
    theta <- array(0, dim = J)
    for (j in 1 : J) {
      for (k in 1:K) {#    
        X <- max(NL[i,j,k,])# all cases are being used to get max NL rating
        theta[j] <- theta[j] + sum( (X < LL[j,,]) * Wk ) + 0.5 * sum( (X == LL[j,,]) * Wk)            
      }  
      theta[j] <- theta[j] / K / K2#
    }
  } 
 
  return (theta)  
}
