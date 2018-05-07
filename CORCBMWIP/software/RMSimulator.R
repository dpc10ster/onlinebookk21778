RMSimulator <- function(I, J, K1, K2, mu, tau, varComp, isBinned, DesiredNumBins){
  stdC <- sqrt(varComp$varC)
  stdTC <- sqrt(varComp$varTC)
  stdRC <- sqrt(varComp$varRC)
  stdEps <- sqrt(varComp$varEps)
  stdR <- sqrt(varComp$varR)
  stdTR <- sqrt(varComp$varTR)
  K <- K1 + K2
  
  Rjt <- array(rnorm(J * 2, sd = stdR), c(J, 2))
  Ckt <- rnorm(K, sd = stdC)
  TRijt <- array(rnorm(I * J * 2, sd = stdTR), c(I, J, 2))
  TCikt <- array(rnorm(I * K, sd = stdTC), c(I, K))
  RCjkt <- array(rnorm(J * K, sd = stdRC), c(J, K))
  EPSijkt <- array(rnorm(I * J * K, sd = stdEps), c(I, J, K))
  
  z1b <- array(dim = c(I, J, K1))
  z2b <- array(dim = c(I, J, K2))
  tempNL <- array(dim = c(I, K1))
  tempLL <- array(dim = c(I, K2))
  for (j in 1:J){
    for (i in 1:I){
      tempNL[i,] <- mu[1]+tau[i, 1] + Rjt[j, 1] + Ckt[1:K1] + 
        TRijt[i, j, 1] + TCikt[i, 1:K1] + RCjkt[j, 1:K1] + EPSijkt[i, j, 1:K1]
      tempLL[i,] <- mu[2]+tau[i, 2] + Rjt[j, 2] + Ckt[(K1 + 1):K] + 
        TRijt[i, j, 2] + TCikt[i, (K1 + 1):K] + RCjkt[j, (K1 + 1):K] + EPSijkt[i, j, (K1 + 1):K]
    }
    if (isBinned){ #  for each reader, must apply same binning rule to both modalities
      x1 <- c(tempNL[1,],tempNL[2,])
      x2 <- c(tempLL[1,],tempLL[2,])
      zb <- ToIntegerRatings(x1, x2, DesiredNumBins)
      z1b[1, j, ] <- zb$f[1:K1]
      z1b[2, j, ] <- zb$f[(K1 + 1):K]
      z2b[1, j, ] <- zb$t[1:K2]
      z2b[2, j, ] <- zb$t[(K2 + 1):K]
    } else {
      z1b[, j, ] <- tempNL
      z2b[, j, ] <- tempLL
    }
  }
  
  return(list(
    z1b = z1b,
    z2b = z2b
  ))
}