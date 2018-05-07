smAuc <- function (mu, lambdaP, nuP, nLesDistr)
{
  maxFPF <- xROC(-Inf, lambdaP)
  AUC <- integrate(ROCInt, 0, maxFPF, mu = mu, lambdaP = lambdaP, nuP = nuP, pmfLesionDistribution = nLesDistr)$value
  AUC <- AUC + (1 + yROC(-Inf, mu, lambdaP, nuP, nLesDistr)) * (1 - maxFPF) / 2
  return (AUC)
}

ZetaOpt <- function(zeta, lambdaP, fpf){
  prob <- xROC(zeta, lambdaP)
  return(prob - fpf)
}

DiffAucSmMinusAucCbm <- function(mu, lambdaP, nuP, nLesDistr, AUCCbm){
  AUCSm <- smAuc (mu, lambdaP, nuP, nLesDistr)
  return ((AUCSm - AUCCbm))
}

ROCInt <- function(FPF, mu, lambdaP, nuP, pmfLesionDistribution){
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp <= 0] <- .Machine$double.xmin
  zeta <- qnorm(tmp)
  if (is.matrix(pmfLesionDistribution)){
    fl <- pmfLesionDistribution[, 2] / sum(pmfLesionDistribution[, 2])
    TPF <- 0
    for (i in 1:nrow(pmfLesionDistribution)){
      TPF <- TPF + fl[i] * (1 - (1 - nuP/2 + nuP/2  *erf( (zeta - mu) / sqrt(2) ))^pmfLesionDistribution[i, 1] * exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
      #TPF <- TPF + (1 - (1 - nuP/2 + nuP/2  *erf( (zeta - mu) / sqrt(2) ))^pmfLesionDistribution[i, 1] * exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
    }
  }else{
    TPF <- (1 - (1 - nuP/2 + nuP/2  *erf( (zeta - mu) / sqrt(2) ))^pmfLesionDistribution * exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
  }
  return (TPF)
}

