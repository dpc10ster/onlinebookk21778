FindParamFixAuc <- function(mu, lambda, nu, lesionDistribution, RocAuc){
  if (missing(mu) || is.null(mu) || is.na(mu)){
    if (AucsRsmRocDiff(0.01, lambda, nu, lesionDistribution, RocAuc) * AucsRsmRocDiff(20, lambda, nu, lesionDistribution, RocAuc) > 0){
      return(NA)
    }
    return(uniroot(AucsRsmRocDiff, c(0.01, 20), lambda = lambda, nu = nu, lesionDistribution = lesionDistribution, RocAuc = RocAuc)$root)
  }
  
  if (missing(lambda) || is.null(lambda) || is.na(lambda)){
    if (AucsRsmRocDiff(mu, 0.01, nu, lesionDistribution, RocAuc) * AucsRsmRocDiff(mu, 20, nu, lesionDistribution, RocAuc) > 0){
      return(NA)
    }
    return(uniroot(AucsRsmRocDiff, c(0.01, 20), mu = mu, nu = nu, lesionDistribution = lesionDistribution, RocAuc = RocAuc)$root)
  }
  
  if (missing(nu) || is.null(nu) || is.na(nu)){
    if (AucsRsmRocDiff(mu, lambda, 0.01, lesionDistribution, RocAuc) * AucsRsmRocDiff(mu, lambda, 1, lesionDistribution, RocAuc) > 0){
      return(NA)
    }
    return(uniroot(AucsRsmRocDiff, c(0.01, 1), mu = mu, lambda = lambda, lesionDistribution = lesionDistribution, RocAuc = RocAuc)$root)
  }
}

AucsRsmRocDiff <- function(mu, lambda, nu, lesionDistribution, RocAuc){
  if (!all(c(length(mu) == length(lambda), length(mu) == length(nu))))
    stop("Parameters have different lengths.")
  
  if (missing(lesionDistribution) || length(lesionDistribution) == 2){
    lesionDistribution <- c(1, 1)
    dim(lesionDistribution) <- c(1, 2)
  }
  
  plotStep <- 0.005
  zeta <- seq(from = -20, to = 20, by = plotStep)
  
  aucROC <- rep(NA, length(mu))
  lambdaP <- lambda
  nuP <- nu
  for (i in 1:length(mu)){
    if (nu[i] < 0 ) stop("nu must be non-negative")
    
    lambdaP[i] <- lambda[i] / mu[i]
    if (abs(nu[i] * mu[i]) <= 1e-6 ) nuP[i] <- 1e-6 else nuP[i] <- (1-exp(-nu[i] * mu[i]))
    maxFPF <- xROC(-20, lambdaP[i])
    maxTPF <- yROC(-20, mu[i], lambdaP[i], nuP[i], lesionDistribution)
    AUC <- integrate(intROC, 0, maxFPF, mu = mu[i], lambdaP = lambdaP[i], nuP = nuP[i], lesionDistribution = lesionDistribution)$value
    aucROC[i] <- AUC + (1 + maxTPF) * (1 - maxFPF) / 2
  }
  return(aucROC - RocAuc)
}