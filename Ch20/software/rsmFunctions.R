
erf <- function(x){
  return (2 * pnorm(sqrt(2) * x) - 1)
}

xROC <- function(zeta, lambdaP){
  # returns FPF, the abscissa of ROC curve
  return( 1 - exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
}

yROC <- function(zeta, mu, lambdaP, nuP, pmfLesionDistribution){
  # returns TPF, the ordinate of ROC curve
  fl <- pmfLesionDistribution[, 2] / sum(pmfLesionDistribution[, 2])
  TPF <- 0
  for (i in 1:nrow(pmfLesionDistribution)){
    TPF <- TPF + fl[i] * (1 - (1 - nuP/2 + nuP/2  *erf( (zeta - mu) / sqrt(2) ))^pmfLesionDistribution[i, 1] * exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2))))
  }
  return (TPF)
}

intROC <- function(FPF, mu, lambdaP, nuP, pmfLesionDistribution){
  # returns TPF, the ordinate of ROC curve; takes FPF as the variable. 
  # AUC is calculated by integrating this function in terms of FPF
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp < 0] <- pnorm(-20)
  zeta <- qnorm(tmp)
  TPF <- yROC(zeta, mu, lambdaP, nuP, pmfLesionDistribution)
  return (TPF)
}

xFROC <- function(zeta, lambdaP){
  # returns NLF, the abscissa of FROC curve
  NLF <- lambdaP * (1 - pnorm(zeta))
  return(NLF)
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

yWAFROC <- function(zeta, mu, nuP, pmfLesionDistribution, lesionWeights){
  # returns wLLFL, the ordinate of wAFROC curve
  fl <- pmfLesionDistribution[, 2] / sum(pmfLesionDistribution[, 2])
  wLLF <- 0
  for (L in 1:nrow(pmfLesionDistribution)){
    nLesion <- pmfLesionDistribution[L, 1] 
    # nLesion is the first element in the row L of pmfLesionDistribution, 
    # which is the number of lesions for this lesion weights distributions condition
    wLLFTmp <- 0
    for (l in 1:nLesion){
      # l is the number of sucesses with number of lesions nLesion
      wLLFTmp <- wLLFTmp + sum(lesionWeights[L, 1:l]) * dbinom(l, nLesion, nuP) * (1 - pnorm(zeta - mu))
      
    }
    wLLF <- wLLF + fl[L] * wLLFTmp
  }
  return(wLLF)
}

intWAFROC <- function(FPF, mu, lambdaP, nuP, pmfLesionDistribution, lesionWeights){
  # returns wLLF, the ordinate of AFROC curve; takes FPF as the variable. 
  # AUC is calculated by integrating this function in terms of FPF
  tmp <- 1 / lambdaP * log(1 - FPF) + 1
  tmp[tmp < 0] <- pnorm(-20)
  zeta <- qnorm(tmp)
  wLLF <- sapply(zeta, yWAFROC, mu = mu, nuP = nuP, pmfLesionDistribution, lesionWeights)
  return(wLLF)
}
