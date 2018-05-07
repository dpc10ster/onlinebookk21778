source("Transforms.R")
source("LL.R")
source("RocOperatingPointsFromRatingsTable.R")
source("VarianceAz.R");
source("ChisqrGoodnessOfFit.R")
source("RocOperatingPointsFromRatingsTable.R")

RocfitR  <- function(K)
{
  
  K <- FixRocDataTable(K)
  if (K[1] == -1) return(-1)
  
  K1 <- K[1,] # this is the observed data!
  K2 <- K[2,] # this is the observed data!
  
  # initial estimates of a and b parameters
  ret <- RocOperatingPointsFromRatingsTable (K1, K2)
  FPF <- ret$FPF; TPF <- ret$TPF; ph_inv_FPF <- qnorm(FPF); ph_inv_TPF <- qnorm(TPF)
  if (length(ph_inv_FPF) == 1) { # if only one point, assume one paramter fit
    return(-1)
    #         a <- (ph_inv_TPF-ph_inv_FPF)
    #         b <- 1
  } else {
    fit <- lm(ph_inv_TPF~ph_inv_FPF) # straight line fit method of estimating a and b  
    a <- fit$coefficients[[1]] # these is the initial estimate of a 
    b <- fit$coefficients[[2]] # these is the initial estimate of b
    if (a > a_max) return(-1)
    if (b > b_max) return(-1)  
    if (a < a_min) return(-1)
    if (b < b_min) return(-1)      
  }
  # thresholds can be estimated by by applying inverse function to Eqn. xx and solving to zeta
  zeta_initial_fpf <- -ph_inv_FPF # see Eqn. xx
  zeta_initial_tpf <- (a - ph_inv_TPF)/b # see Eqn. xx
  zeta_initial <- (zeta_initial_fpf + zeta_initial_tpf)/2 # average the two estimates
  zeta_initial <- rev(zeta_initial) # apply reverse order to correct the ordering of the cutoffs
  zeta_initial_guess <- seq(-b, a + 1, length.out = length(K1)-1) # to test stability of alg. to guess choice
  
  param_initial <- c(a, b, zeta_initial) 
  #param_initial <- c(1, 1, zeta_initial_guess) # to test stability of alg. to other choices
  
  param_initial_prime <- ThetaPrime(param_initial)# use this method to test variation of -LL with parameters
  LLval_initial <- LL(param_initial_prime, K1, K2) # use this method to test variation of -LL with parameters
  
  ret_nlm <- nlm(LL, param_initial_prime, K1 = K1, K2 = K2, stepmax = 1) # this does the actual minimization of -LL
  parameters_final_nlm <- Theta(ret_nlm$estimate)
  
  #   if (length(ph_inv_FPF) != 1) { 
  options(warn=-1)
  hess <- hessian(LL_theta, parameters_final_nlm, method="Richardson", K1 = K1, K2 = K2)
  options(warn=0)
  if (any(is.nan(hess))) return(-1)
  Cov <- solve(hess)
  StdAz <- sqrt(VarianceAz (a, b,Cov))
  #   } else {
  #     Cov <- 0    
  #     StdAz <- 0    
  #   }
  Az <- pnorm(a/sqrt(1+b^2))  
  
  a <- parameters_final_nlm[1];b <- parameters_final_nlm[2]
  zeta <- parameters_final_nlm[3:length(parameters_final_nlm)]
  
  return (list( 
    Az = Az,
    StdAz = StdAz,
    a = a,
    b = b,
    zeta  = zeta,
    Cov = Cov))
}


FixRocDataTable  <- function (RocDataTable) {
  R  <- length(RocDataTable[1,])
  while(1) {
    if (!any(RocDataTable == 0)) return (RocDataTable) else {
      for (i in 1:2) {
        if (!any(RocDataTable[i,] == 0)) next else {
          if (R - length(which(RocDataTable[i,] == 0)) <= 2) return(-1)#this needs some more thought
          newBins  <- R - length(which(RocDataTable[i,] == 0))
          zeroColumn <- which(RocDataTable[i,] == 0)[1]
          if (zeroColumn == R) {
            RocDataTable[1,R-1]  <- sum(RocDataTable[1,(R-1):R]) 
            RocDataTable[2,R-1]  <- sum(RocDataTable[2,(R-1):R])
            RocDataTable  <- RocDataTable[,-zeroColumn]
          } else if (zeroColumn == 1) {
            RocDataTable[1,2]  <- sum(RocDataTable[1,1:2]) 
            RocDataTable[2,2]  <- sum(RocDataTable[2,1:2])
            RocDataTable  <- RocDataTable[,-zeroColumn]
          } else {
            RocDataTable[1,zeroColumn-1]  <- sum(RocDataTable[1,(zeroColumn-1):zeroColumn]) 
            RocDataTable[2,zeroColumn-1]  <- sum(RocDataTable[2,(zeroColumn-1):zeroColumn])
            RocDataTable  <- RocDataTable[,-zeroColumn]         
          }
          R  <- R - 1
        }
      }
    }
    if (!any(RocDataTable == 0)) return (RocDataTable) 
  }
}

