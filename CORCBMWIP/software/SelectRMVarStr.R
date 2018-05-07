# limited to two modalities
SelectRMVarStr <- function(I,AxOrgBar,RMVarStr,nhCondition)
{
  mu <- array(0, dim = 2)
  tau <- array(0, dim = c(2,2))
  Auci <- AxOrgBar
  separations <- sqrt(2)*qnorm(Auci)
  mu[2] <- mean(separations)
  tau[2, 1] <- sqrt(2)*qnorm(Auci[1]) - mu[2]
  tau[2, 2] <- sqrt(2)*qnorm(Auci[2]) - mu[2]
  
  if (RMVarStr == "HL"){
    varComp <- RoeMetzVarStr(1)
  } else if (RMVarStr == "LL"){
    varComp <- RoeMetzVarStr(2)
  } else if (RMVarStr == "HH"){
    if (mu[2] <= sqrt(2)*qnorm(0.702)){
      varComp <- RoeMetzVarStr(3)
    }else if (mu[2] > sqrt(2)*qnorm(0.702) && mu[2] < sqrt(2)*qnorm(0.962)){
      varComp <- RoeMetzVarStr(4)
    }else if (mu[2] >= sqrt(2)*qnorm(0.962)){
      varComp <- RoeMetzVarStr(5)
    }
  } else if (RMVarStr == "LH"){
    if (mu[2] <= sqrt(2)*qnorm(0.702)){
      varComp <- RoeMetzVarStr(6)
    }else if (mu[2] > sqrt(2)*qnorm(0.702) && mu[2] < sqrt(2)*qnorm(0.962)){
      varComp <- RoeMetzVarStr(7)
    }else if (mu[2] >= sqrt(2)*qnorm(0.962)){
      varComp <- RoeMetzVarStr(8)
    }
  } else if (RMVarStr == "EST"){
    varComp <- RmVarEstimate(orgData)
  }
  if (nhCondition == TRUE) {
    tau[2,] <- 0
  }
  
  return(list(
    varComp = varComp,
    mu = mu,
    tau = tau
  ))
}