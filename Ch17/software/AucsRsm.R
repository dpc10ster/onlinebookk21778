AucsRsm <- function(mu, lambda, nu, lesionDistribution){
  if (!all(c(length(mu) == length(lambda), length(mu) == length(nu))))
    stop("Parameters have different lengths.")
  
  if (missing(lesionDistribution)){
    lesionDistribution <- c(1, 1)
    dim(lesionDistribution) <- c(1, 2)
  }
  
  plotStep <- 0.005
  zeta <- seq(from = -20, to = 20, by = plotStep)
  
  aucROC <- rep(NA, length(mu))
  aucAFROC <- aucROC
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
    
    maxLLF <- yFROC(-20, mu[i], nuP[i])
    AUC <- integrate(intAFROC, 0, maxFPF, mu = mu[i], lambdaP = lambdaP[i], nuP = nuP[i])$value
    aucAFROC[i] <- AUC + (1 + maxLLF) * (1 - maxFPF) / 2
  }
  return(list(
    aucROC = aucROC,
    aucAFROC = aucAFROC
  ))
}

erf <- function(x){
  return (2 * pnorm(sqrt(2) * x) - 1)
}

xROC <- function(zeta, lambdaP){
  # returns FPF, the abscissa of ROC curve
  return( 1 - exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
}

yROC <- function(zeta, mu, lambdaP, nuP, lesionDistribution){
  # returns TPF, the ordinate of ROC curve
  fl <- lesionDistribution[, 2] / sum(lesionDistribution[, 2])
  TPF <- 0
  for (i in 1:nrow(lesionDistribution)){
    TPF <- TPF + fl[i] * (1 - (1 - nuP/2 + nuP/2  *erf( (zeta - mu) / sqrt(2) ))^lesionDistribution[i, 1] * exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
  }
  return (TPF)
}

intROC <- function(FPF, mu, lambdaP, nuP, lesionDistribution){
  # returns TPF, the ordinate of ROC curve; takes FPF as the variable. 
  # AUC is calculated by integrating this function in terms of FPF
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp < 0] <- pnorm(-20)
  zeta <- qnorm(tmp)
  TPF <- yROC(zeta, mu, lambdaP, nuP, lesionDistribution)
  return (TPF)
}

yFROC <- function(zeta, mu, nuP){
  # returns LLF, the ordinate of FROC, AFROC curve
  LLF <- nuP * (1 - pnorm(zeta - mu))
  return(LLF)
}

intAFROC <- function(FPF, mu, lambdaP, nuP){
  # returns LLF, the ordinate of AFROC curve; takes FPF as the variable. 
  # AUC is calculated by integrating this function in terms of FPF
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp < 0] <- pnorm(-20)
  zeta <- qnorm(tmp)
  LLF <- yFROC(zeta, mu, nuP)
  return(LLF)
}