ForwardTransform <- function(parametersArray, upperLambda) {
  maxLambdaP <- RJafrocEnv$maxLambdaP
  minLambdaP <- RJafrocEnv$minLambdaP
  maxNuP <- RJafrocEnv$maxNuP
  minNuP <- RJafrocEnv$minNuP
  if (missing(upperLambda)){
    maxLambdaP <- RJafrocEnv$maxLambdaP
  }else{
    maxLambdaP <- upperLambda
  }  
  minLambdaP <- RJafrocEnv$minLambdaP
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
  maxNuP <- RJafrocEnv$maxNuP
  minNuP <- RJafrocEnv$minNuP
  if (missing(upperLambda)){
    maxLambdaP <- RJafrocEnv$maxLambdaP
  }else{
    maxLambdaP <- upperLambda
  }
  minLambdaP <- RJafrocEnv$minLambdaP
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
  maxNuP <- RJafrocEnv$maxNuP
  minNuP <- RJafrocEnv$minNuP
  nuP <- 1 - exp(-nu * abs(mu))
  if (nuP > (1 - minNuP)) 
    nuP <- 1 - minNuP
  if (nuP < minNuP) 
    nuP <- minNuP
  return(nuP)
} 

ForwardLambda <- function(lambdaP, lambdPLower){
  maxLambdaP <- RJafrocEnv$maxLambdaP
  minLambdaP <- RJafrocEnv$minLambdaP
  return(log(-log((lambdaP - minLambdaP)/(maxLambdaP - minLambdaP))))
}

InverseLamdba <- function(lambdaP, lambdPLower){
  maxLambdaP <- RJafrocEnv$maxLambdaP
  minLambdaP <- RJafrocEnv$minLambdaP
  return(minLambdaP + (maxLambdaP - minLambdaP) * exp(-exp(lambdaP)))
}

ForwardNu <- function(nuP, nuPLower, nuPUpper){
  maxNuP <- RJafrocEnv$maxNuP
  minNuP <- RJafrocEnv$minNuP
  return(log(-log((nuP - minNuP)/(maxNuP - minNuP))))
}

InverseNu <- function(nuP, nuPLower, nuPUpper){
  maxNuP <- RJafrocEnv$maxNuP
  minNuP <- RJafrocEnv$minNuP
  return(minNuP + (maxNuP - minNuP) * exp(-exp(nuP)))
}

ForwardZetas <- function(zetas){
  zetasFwd <- zetas
  if (length(zetas) > 1){
    zetasFwd[2:length(zetasFwd)] <- log(zetasFwd[2:length(zetasFwd)] - zetasFwd[1:(length(zetasFwd) - 1)])
  }
  # zetasFwd[1] <- ForwardValue(zetasFwd[1], -4, 3)
  # for (z in 2:length(zetas))
  # zetasFwd[z] <- ForwardValue(zetasFwd[z], -7, maxZeta + mu - zetas)
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