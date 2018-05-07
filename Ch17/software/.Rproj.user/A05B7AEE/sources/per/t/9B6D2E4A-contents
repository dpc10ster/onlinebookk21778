# MainRsmSwetsObservations.R
rm(list = ls())
library(RJafroc)
source("rsmPdfMeansAndStddevs.R")
source("FindParamFixAuc.R")

logseq <- function( d1, d2, n) {
  logf <- log(d2/d1)/(n-1)
  return (exp(seq(log(d1), log(d2), logf)))
}

Lmax <- 1;K2 <- 700;Lk2 <- floor(runif(K2, 1, Lmax + 1))
nLesPerCase <- unique(Lk2);lesionDist <- array(dim = c(length(nLesPerCase), 2))
for (i in nLesPerCase) lesionDist[i, ] <- c(i, sum(Lk2 == i)/K2)

RsmRocAuc <- 0.8 # other parameters adjusted to attain this value; 0.7 or 0.8
cat("RsmRocAuc constraint = ", RsmRocAuc, "\n")

# Part I
lambda <- 2;cat("\nVary mu and nu only, ", "lambda = ", lambda, "\n")
muArr <- logseq(2,5,10)
for (i in 1:length(muArr)) {
  mu <- muArr[i];nu <- NA # intrinsic parameters
  retParms <- FindParamFixAuc(mu, lambda, nu, lesionDist, RsmRocAuc)
  if (!is.na(retParms)) nu <- retParms else next
  ret <- rsmPdfMeansAndStddevs(mu, lambda, nu, lesionDist)
  meanN <- ret$meanN;meanD <- ret$meanD;stdDevN <- ret$stdErrN;stdDevD <- ret$stdErrD

  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, type = "ROC", 
    lesionDistribution = lesionDist, legendPosition  = "none")
  cat("mu = ", mu,", lambda = ", lambda,
      ", nu = ", nu, ", AUC = ", ret1$aucROC,
      ", bParm = ", stdDevN/stdDevD, 
      ", dmu/dsigma = ",  (meanD - meanN)/(stdDevD - stdDevN), "\n")
  next
}

# Part II
nu <- 1;cat("\nVary mu and lambda only, ", "nu = ", nu, "\n")
lambdaArr <- logseq(1, 5, 10)
for (i in 1:length(lambdaArr)) {
  lambda <- lambdaArr[i]; mu <- NA; # intrinsic parameters

  retParms <- FindParamFixAuc(mu, lambda, nu, lesionDist, RsmRocAuc)
  if (!is.na(retParms)) mu <- retParms else next
  ret <- rsmPdfMeansAndStddevs(mu, lambda, nu, lesionDist)
  meanN <- ret$meanN;meanD <- ret$meanD
  stdDevN <- ret$stdErrN;stdDevD <- ret$stdErrD

  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, type = "ROC", 
    lesionDistribution = lesionDist, legendPosition  = "none")
  cat("mu = ", mu,", lambda = ", lambda,
      ", nu = ", nu, ", AUC = ", ret1$aucROC,
      ", bParm = ", stdDevN/stdDevD, 
      ", dmu/dsigma = ",  (meanD - meanN)/(stdDevD - stdDevN), "\n")
  next
}

# Part III
mu <- 2;cat("\nVary lambda and nu only,", "mu = ", mu, "\n")
lambdaArr <- logseq(0.1, 5, 10)
for (i in 1:length(lambdaArr)) {
  lambda <- lambdaArr[i]; nu <- NA;  # intrinsic parameters
  retParms <- FindParamFixAuc(mu, lambda, nu, lesionDist, RsmRocAuc)
  if (!is.na(retParms)) nu <- retParms else next
  ret <- rsmPdfMeansAndStddevs(mu, lambda, nu, lesionDist)
  meanN <- ret$meanN;meanD <- ret$meanD
  stdDevN <- ret$stdErrN;stdDevD <- ret$stdErrD
  
  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, type = "ROC", 
    lesionDistribution = lesionDist, legendPosition  = "none")
  cat("mu = ", mu,", lambda = ", lambda,
      ", nu = ", nu, ", AUC = ", ret1$aucROC,
      ", bParm = ", stdDevN/stdDevD, 
      ", dmu/dsigma = ",  (meanD - meanN)/(stdDevD - stdDevN), "\n")
  next
}
