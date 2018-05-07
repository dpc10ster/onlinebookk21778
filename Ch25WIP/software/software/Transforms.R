ForwardTransform <- function(parametersArray, upperLambda) {
  if (!missing(upperLambda)){
    maxLambdaP <- upperLambda
  }
  if (length(parametersArray) == 1) {
    parametersTransformed <- log(-log((parametersArray - minLambdaP)/(maxLambdaP - minLambdaP)))
  } else {
    parametersTransformed <- parametersArray
    parametersTransformed[1] <- log(-log((parametersArray[1] - minLambdaP)/(maxLambdaP - minLambdaP)))
    parametersTransformed[2] <- log(-log((parametersArray[2] - minNuP)/(maxNuP - minNuP)))
    if (length(parametersArray) > 3) {
      parametersTransformed[4:length(parametersArray)] <- log(parametersArray[4:length(parametersArray)] - parametersArray[3:(length(parametersArray) - 1)])
      parametersTransformed[4:length(parametersArray)] <- apply(rbind(parametersTransformed[4:length(parametersArray)], -2000), 2, max)
    }
  }
  return(parametersTransformed)
}

InverseTransform <- function(parametersTransformed, upperLambda) {
  if (!missing(upperLambda)){
    maxLambdaP <- upperLambda
  }
  if (length(parametersTransformed) == 1) {
    parametersArray <- minLambdaP + (maxLambdaP - minLambdaP) * exp(-exp(parametersTransformed))
  } else {
    parametersArray <- parametersTransformed
    parametersArray[1] <- minLambdaP + (maxLambdaP - minLambdaP) * exp(-exp(parametersTransformed[1]))
    parametersArray[2] <- minNuP + (maxNuP - minNuP) * exp(-exp(parametersTransformed[2]))
    if (length(parametersTransformed) > 3) {
      for (i in 4:length(parametersTransformed)) parametersArray[i] <- exp(parametersTransformed[i]) + parametersArray[i - 1]
    }
  }
  return(parametersArray)
}

NuToNuP <- function(nu, mu) {
  nuP <- 1 - exp(-nu * abs(mu))
  if (nuP > (1 - 1e-07)) 
    nuP <- 1 - 1e-07
  if (nuP < 1e-07) 
    nuP <- 1e-07
  return(nuP)
} 

ForwardLambda <- function(lambdaP, lambdPLower){
  minLambdaP <- lambdPLower
  return(log(-log((lambdaP - minLambdaP)/(maxLambdaP - minLambdaP))))
}

InverseLamdba <- function(lambdaP, lambdPLower){
  minLambdaP <- lambdPLower
  return(minLambdaP + (maxLambdaP - minLambdaP) * exp(-exp(lambdaP)))
}

ForwardNu <- function(nuP, nuPLower, nuPUpper){
  maxNuP <- nuPUpper
  minNuP <- nuPLower
  return(log(-log((nuP - minNuP)/(maxNuP - minNuP))))
}

InverseNu <- function(nuP, nuPLower, nuPUpper){
  maxNuP <- nuPUpper
  minNuP <- nuPLower
  return(minNuP + (maxNuP - minNuP) * exp(-exp(nuP)))
}

ForwardZetas <- function(zetas){
  zetasFwd <- zetas
  if (length(zetas) > 1){
    zetasFwd[2:length(zetasFwd)] <- log(zetasFwd[2:length(zetasFwd)] - zetasFwd[1:(length(zetasFwd) - 1)])
  }
  # zetasFwd[1] <- ForwardValue(zetasFwd[1], -4, 3)
  # zetasFwd[-1] <- ForwardValue(zetasFwd[-1], -7, 1)
  return(zetasFwd)
}

InverseZetas <- function(zetasFwd){
  # zetasFwd[1] <- InverseValue(zetasFwd[1], -4, 3)
  # zetasFwd[-1] <- InverseValue(zetasFwd[-1], -7, 1)
  zetas <- zetasFwd
  if (length(zetasFwd) > 1) {
    for (i in 2:length(zetasFwd)) zetas[i] <- exp(zetasFwd[i]) + zetas[i - 1]
  }
  return(zetas)
}

ForwardValue <- function(value, valueLower, valueUpper){
  maxValue <- valueUpper
  minValue <- valueLower
  return(log(-log((value - minValue)/(maxValue - minValue))))
}

InverseValue <- function(valueFwd, valueLower, valueUpper){
  maxValue <- valueUpper
  minValue <- valueLower
  return(minValue + (maxValue - minValue) * (exp(-exp(valueFwd))))
}