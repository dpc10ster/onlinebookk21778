source("AucsRsm.R")
rsmPdfMeansAndStddevs <- function(mu, lambda, nu, lesionDistribution){
  lambdaP <- lambda / mu
  if (abs(nu * mu) <= 1e-6 ) nuP <- 1e-6 else nuP <- (1-exp(-nu * mu))
  FPFmax <- xROC(-Inf, lambdaP)
  meanN <- integrate(meanNInt, -Inf, Inf, lambdaP = lambdaP, FPFmax = FPFmax)$value
  varN <- integrate(varNInt, -Inf, Inf, lambdaP = lambdaP, FPFmax = FPFmax)$value - meanN^2 # https://en.wikipedia.org/wiki/Variance
  stdErrN <- sqrt(varN)
  TPFmax <- yROC(-Inf, mu, lambdaP, nuP, lesionDistribution)
  meanD <- integrate(meanDInt, -Inf, Inf, mu = mu, lambdaP = lambdaP, nuP = nuP, lesionDistribution = lesionDistribution, TPFmax = TPFmax)$value
  varD <- integrate(varDInt, -Inf, Inf, mu = mu, lambdaP = lambdaP, nuP = nuP, lesionDistribution = lesionDistribution, TPFmax = TPFmax)$value - meanD^2 # https://en.wikipedia.org/wiki/Variance
  stdErrD <- sqrt(varD)
  return(list(
    meanN = meanN,
    stdErrN = stdErrN,
    meanD = meanD,
    stdErrD = stdErrD
  ))
}

pdfN <- function(zeta, lambdaP){
  # verified using pdfN(1, 1) and (xROC(1, 1) - xROC(1 + 1e-8, 1)) / 1e-8
  return(lambdaP * exp(-zeta^2/2) * exp(-lambdaP/2 * (1 - erf(zeta/sqrt(2)))) / sqrt(2 * pi))
}

meanNInt <- function(zeta, lambdaP, FPFmax){
  return(pdfN(zeta, lambdaP) / FPFmax * zeta)
}

varNInt <- function(zeta, lambdaP, FPFmax){
  return(pdfN(zeta, lambdaP) / FPFmax * zeta^2)
}

pdfD <- function(zeta, mu, lambdaP, nuP, lesionDistribution){
  # verified using pdfD(1, 1, 1, 0.5, rbind(c(1, 50), c(2, 30))) and 
  # (yROC(1, 1, 1, 0.5, rbind(c(1, 50), c(2, 30))) - yROC(1 + 1e-8, 1, 1, 0.5, rbind(c(1, 50), c(2, 30)))) / 1e-8
  fl <- lesionDistribution[, 2] / sum(lesionDistribution[, 2])
  pdfTemp <- 0
  for (i in 1:nrow(lesionDistribution)){
    pdfTemp <- pdfTemp + fl[i] * (((1 - nuP/2 + nuP/2  *erf( (zeta - mu) / sqrt(2) ))^(lesionDistribution[i, 1] - 1) * exp( (-lambdaP / 2) + 0.5 * lambdaP * erf(zeta / sqrt(2)))) *
                                    (lesionDistribution[i, 1] * nuP/sqrt(2 * pi) * exp(-(zeta - mu)^2/2) + (1 - nuP/2 + nuP/2  *erf( (zeta - mu) / sqrt(2) )) * lambdaP/sqrt(2*pi) * exp(-zeta^2/2)))
  }
  return (pdfTemp)
}

meanDInt <- function(zeta, mu, lambdaP, nuP, lesionDistribution, TPFmax){
  return(pdfD(zeta, mu, lambdaP, nuP, lesionDistribution) / TPFmax * zeta)
}

varDInt <- function(zeta, mu, lambdaP, nuP, lesionDistribution, TPFmax){
  return(pdfD(zeta, mu, lambdaP, nuP, lesionDistribution) / TPFmax * zeta^2)
}