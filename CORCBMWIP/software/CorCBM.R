CorCBM <- function(FPCounts, TPCounts, aucX, aucY, K1, K2, rhoIniNor, rhoIniAbn){
  if (pnorm(muMax / sqrt(2)) <= max(c(aucX, aucY))){
    muMax <<- qnorm(max(c(aucX, aucY))) * sqrt(2) + 0.5
    zetaMax <<- muMax + 2
  }
  
  zeroBinsX <- intersect(which(rowSums(FPCounts) == 0),  which(rowSums(TPCounts) == 0))
  if (length(zeroBinsX) > 0){
    FPCounts <- FPCounts[-zeroBinsX, ]
    TPCounts <- TPCounts[-zeroBinsX, ]
  }
  iniParamX <- IniParam(rowSums(FPCounts), rowSums(TPCounts), aucX, K1, K2)
  muIniX <- iniParamX$muIni
  zetaIniX <- iniParamX$zetaIni
  alphaIniX <- iniParamX$alphaIni
  
  zeroBinsY <- intersect(which(colSums(FPCounts) == 0),  which(colSums(TPCounts) == 0))
  if (length(zeroBinsY) > 0){
    FPCounts <- FPCounts[ , -zeroBinsY]
    TPCounts <- TPCounts[ , -zeroBinsY]
  }
  iniParamY <- IniParam(colSums(FPCounts), colSums(TPCounts), aucY, K1, K2)
  muIniY <- iniParamY$muIni
  zetaIniY <- iniParamY$zetaIni
  alphaIniY <- iniParamY$alphaIni
  
  muIniXFwd <- ForwardValue(muIniX, muMin, muMax)
  muIniYFwd <- ForwardValue(muIniY, muMin, muMax)
  alphaIniXFwd <- ForwardValue(alphaIniX, alphaMin, alphaMax)
  alphaIniYFwd <- ForwardValue(alphaIniY, alphaMin, alphaMax)
  rhoIniNorFwd <- ForwardValue(rhoIniNor, rhoMin, rhoMax)
  rhoIniAbn2Fwd <-  ForwardValue(rhoIniAbn, rhoMin, rhoMax)
  zetaIniXFwd <- ForwardZetas(zetaIniX)
  zetaIniYFwd <- ForwardZetas(zetaIniY)
  
  parameters <- c(list(muIniXFwd, muIniYFwd, alphaIniXFwd, alphaIniYFwd, rhoIniNorFwd, rhoIniAbn2Fwd), as.list(zetaIniXFwd), as.list(zetaIniYFwd))
  # parameters <- c(list(muIniXFwd, muIniYFwd, alphaIniXFwd, alphaIniYFwd, ForwardValue(0.1, rhoMin, rhoMax), ForwardValue(0.8, rhoMin, rhoMax)), as.list(zetaIniXFwd), as.list(zetaIniYFwd))
  namesVector <- c("muXFwd", "muYFwd", "alphaXFwd", "alphaYFwd","rhoNorFwd", "rhoAbn2Fwd")
  for (z in 1:length(zetaIniX)){
    namesVector <- c(namesVector, paste0("zetaXFwd", z))
  }
  for (z in 1:length(zetaIniY)){
    namesVector <- c(namesVector, paste0("zetaYFwd", z))
  }
  names(parameters) <- namesVector
  
  CorCBMNLLNew <- AddZetaXFwd(CorCBMNLL, length(zetaIniX))
  CorCBMNLLNew <- AddZetaYFwd(CorCBMNLLNew, length(zetaIniY))
  # ret <- mle2(CorCBMNLLTest, start = parameters, method = "BFGS", data = list(FPCounts = FPCounts, TPCounts = TPCounts))
  ret <- mle2(CorCBMNLLNew, start = parameters, method = "BFGS", data = list(FPCounts = FPCounts, TPCounts = TPCounts))
  
  allCoef <- ret@coef
  muX <- as.numeric(InverseValue(allCoef[1], muMin, muMax))
  muY <- as.numeric(InverseValue(allCoef[2], muMin, muMax))
  alphaX <- as.numeric(InverseValue(allCoef[3], alphaMin, alphaMax))
  alphaY <- as.numeric(InverseValue(allCoef[4], alphaMin, alphaMax))
  rhoNor <- as.numeric(InverseValue(allCoef[5], rhoMin, rhoMax))
  rhoAbn2 <- as.numeric(InverseValue(allCoef[6], rhoMin, rhoMax))
  zetaX <- as.numeric(InverseZetas(allCoef[7:(6 + length(zetaIniX))]))
  zetaY <- as.numeric(InverseZetas(allCoef[(7 + length(zetaIniX)):length(allCoef)]))
  
  param <- c(muX, muY, alphaX, alphaY, rhoNor, rhoAbn2, zetaX, zetaY)
  fixParam <- 2 + which(c(alphaX, alphaY, abs(rhoNor), abs(rhoAbn2)) > 0.99)
  if (length(fixParam) > 0){
    param <- param[-fixParam]
    for (p in fixParam){
      if (p == 3){
        alphaX <- 1
      }else if (p == 4){
        alphaY <- 1
      }else if (p == 5){
        rhoNor <- ifelse(test = rhoNor > 0, yes = 1, no = -1)
      }else if (p == 6){
        rhoAbn2 <- ifelse(test = rhoAbn2 > 0, yes = 1, no = -1)
      }
    }
  }
  
  hess <- hessian(CorCBMLLNoTrns, param, method.args = list(d = 1e-3), FPCounts = FPCounts, TPCounts = TPCounts, fixParam = fixParam)
  
  paramLength <- 6 - length(fixParam)
  covMat <- solve(-hess)[1:paramLength, 1:paramLength]
  vars <- diag(covMat)
  stdErr <- sqrt(vars)
  dMuX <- alphaX * dnorm(muX / sqrt(2)) / sqrt(2)
  dAlphaX <- -1/2 + pnorm(muX / sqrt(2))
  stdAucX <- sqrt(vars[1] * dMuX^2 + vars[3] * dAlphaX^2 + dMuX * dAlphaX * covMat[1, 3])
  
  dMuY <- alphaY * dnorm(muY / sqrt(2)) / sqrt(2)
  dAlphaY <- -1/2 + pnorm(muY / sqrt(2))
  stdAucY <- sqrt(vars[2] * dMuY^2 + vars[4] * dAlphaY^2 + dMuY * dAlphaY * covMat[2, 4])
  
  stdErr <- c(stdErr, stdAucX, stdAucY)
  
  if (length(fixParam) > 0){
    for (p in fixParam){
      stdErr <- append(stdErr, NA, after = p - 1)
      if (p == 3){
        alphaX <- 1
      }else if (p == 4){
        alphaY <- 1
      }else if (p == 5){
        rhoNor <- 1
      }else if (p == 6){
        rhoAbn2 <- 1
      }
    }
    stdErr <- as.numeric(stdErr)
  }
  stdErr[2:3] <- stdErr[3:2] # switch the position of alphaX and muY
  
  
  aucX <- (1 - alphaX) * 0.5 + alphaX * pnorm(muX / sqrt(2))
  aucY <- (1 - alphaY) * 0.5 + alphaY * pnorm(muY / sqrt(2))
  
  aucDiff <- aucX - aucY
  
  varDiff <- 0
  derivs <- c(dMuX, -dMuY, dAlphaX, -dAlphaY)
  for (k in 1:4){
    for (l in 1:4){
      varDiff <- derivs[l] * derivs[k] * covMat[k, l] + varDiff
    }
  }
  stdDiff <- sqrt(varDiff)
  aStat <- abs(aucDiff) / stdDiff
  aPval <- 1 - pnorm(aStat)
  
  return(list(muX = muX, 
              alphaX = alphaX,
              muY = muY, 
              alphaY = alphaY,
              rhoNor = rhoNor, 
              rhoAbn2 = rhoAbn2, 
              zetaX = zetaX, 
              zetaY = zetaY, 
              aucX = aucX,
              aucY = aucY,
              stdErr = stdErr, 
              aStat = aStat,
              aPval = aPval, 
              covMat = covMat))
}

IniParam <- function(fpCounts, tpCounts, auc, K1, K2){
  tpf <- cumsum(rev(tpCounts))/sum(tpCounts)
  fpf <- cumsum(rev(fpCounts))/sum(fpCounts)
  alphaIni <- rev(tpf)[2]
  # alphaIni <- 1 - (1 - rev(tpf)[2]) / (1 - rev(fpf)[2])
  
  muIni <- min(qnorm(min((auc - (1 - alphaIni)/2) / alphaIni, 1)) * sqrt(2), muMax)
  fpCountsCum <- cumsum(fpCounts)[1:(length(fpCounts) - 1)]
  fpCountsCum[fpCountsCum == 0] <- pnorm(-zetaMin)
  tpCountsCum <- cumsum(tpCounts)[1:(length(tpCounts) - 1)]
  tpCountsCum[tpCountsCum == 0] <- pnorm(-zetaMin, mean = muIni)
  
  zetaNorIni <- qnorm(fpCountsCum / K1)
  zetaAbnIni <- sapply(tpCountsCum / K2, IniZetaAbn, alpha = alphaIni, mu = muIni)
  for (z in 1:length(zetaNorIni)){
    zetaNorIni[z] <- min(max(zetaNorIni[z], zetaNorIni[z - 1] + 0.001), zetaMax)
    zetaAbnIni[z] <- min(max(zetaAbnIni[z], zetaAbnIni[z - 1] + 0.001), zetaMax)
  }
  zetaIni <- (zetaNorIni + zetaAbnIni) / 2
  # zetaIni <- zetaNorIni
  
  
  muIni <- min(muIni, muMax - 0.001)
  alphaIni <- min(alphaIni, alphaMax - 0.001)
  # muIniFwd <- ForwardValue(muIni, muMin, muMax)
  # alphaIniFwd <- ForwardValue(alphaIni, alphaMin, alphaMax)
  # zetaIniFwd <- ForwardZetas(zetaIni)
  # parameters <- c(list(muIniFwd, alphaIniFwd), as.list(zetaIniFwd))
  # namesVector <- c("muFwd", "alphaFwd")
  # for (z in 1:length(zetaIni)){
  #   namesVector <- c(namesVector, paste0("zetaFwd", z))
  # }
  # names(parameters) <- namesVector
  # CBMNLLNew <- AddZetaFwd(CBMNLL, length(zetaIni))
  # ret <- mle2(CBMNLLNew, start = parameters, method = "BFGS", data = list(FPCounts = fpCounts, TPCounts = tpCounts))
  # allCoef <- ret@coef
  # muFwd <- allCoef[1]
  # alphaFwd <- allCoef[2]
  # zetaFwd <- allCoef[3:length(allCoef)]
  # 
  # muIni <- InverseValue(muFwd, muMin, muMax)
  # alphaIni <- InverseValue(alphaFwd, alphaMin, alphaMax)
  # zetaIni <- InverseZetas(zetaFwd)
  
  return(list(muIni = muIni,
              zetaIni = zetaIni, 
              alphaIni = alphaIni))
}

IniZetaAbn <- function(alpha, mu, tpCum){
  return(uniroot(TPCumDiff, c(zetaMin, zetaMax), alpha = alpha, mu = mu, tpCum = tpCum)$root)
}

TPCumDiff <- function(zeta, alpha, mu, tpCum){
  return((1 - alpha) * pnorm(zeta) + alpha * pnorm(zeta - mu) - tpCum)
}

CBMNLL <- function(muFwd, alphaFwd){
  allParameters <- names(formals())
  zetaPos <- regexpr("zeta", allParameters)
  zetaFwd <- unlist(mget(allParameters[which(zetaPos == 1)]))
  
  mu <- InverseValue(muFwd, muMin, muMax)
  alpha <- InverseValue(alphaFwd, alphaMin, alphaMax)
  zeta <- InverseZetas(zetaFwd)
  
  zeta <- c(-Inf, zeta, Inf)
  
  LLNor <- 0
  LLAbn <- 0
  nBins <- length(FPCounts)
  for (b in 1:nBins){
    LLNor <- LLNor + FPCounts[b] * log(pnorm(zeta[b + 1]) - pnorm(zeta[b]))
    LLAbn <- LLAbn + TPCounts[b] * log((1 - alpha) * (pnorm(zeta[b + 1]) - pnorm(zeta[b])) + alpha * (pnorm(zeta[b + 1], mean = mu) - pnorm(zeta[b], mean = mu)))
  }
  LL <- LLNor + LLAbn
  if (any(is.infinite(c(LLNor, LLAbn)))){
    Stop <- 1
  }
  
  return(-LL)
}

CorCBMNLL <- function(muXFwd, muYFwd, alphaXFwd, alphaYFwd, rhoNorFwd, rhoAbn2Fwd){
  allParameters <- names(formals())
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaXFwd <- unlist(mget(allParameters[which(zetaXPos == 1)]))
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaYFwd <- unlist(mget(allParameters[which(zetaYPos == 1)]))
  
  muX <- InverseValue(muXFwd, muMin, muMax)
  muY <- InverseValue(muYFwd, muMin, muMax)
  alphaX <- InverseValue(alphaXFwd, alphaMin, alphaMax)
  alphaY <- InverseValue(alphaYFwd, alphaMin, alphaMax)
  rhoNor <- InverseValue(rhoNorFwd, rhoMin, rhoMax)
  rhoAbn1 <- InverseValue(rhoNorFwd, rhoMin, rhoMax)
  rhoAbn2 <- InverseValue(rhoAbn2Fwd, rhoMin, rhoMax)
  rhoAbn3 <- mean(c(rhoAbn1, rhoAbn2))
  rhoAbn4 <- rhoAbn3
  zetaX <- InverseZetas(zetaXFwd)
  zetaY <- InverseZetas(zetaYFwd)
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaAbn1 <- rbind(c(1, rhoAbn1), c(rhoAbn1, 1))
  sigmaAbn2 <- rbind(c(1, rhoAbn2), c(rhoAbn2, 1))
  sigmaAbn3 <- rbind(c(1, rhoAbn3), c(rhoAbn3, 1))
  sigmaAbn4 <- rbind(c(1, rhoAbn4), c(rhoAbn4, 1))
  
  LLNor <- 0
  LLAbn <- 0
  nBinsX <- dim(FPCounts)[1]
  nBinsY <- dim(FPCounts)[2]
  for (bX in 1:nBinsX){
    for (bY in 1:nBinsY){
      LLNor <- LLNor + FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
      LLAbn <- LLAbn + TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
    }
  }
  LL <- LLNor + LLAbn
  
  return(-LL)
}


CorCBMLLNoTrns <- function(param, FPCounts, TPCounts, fixParam){
  if (length(fixParam) > 0){
    for (p in fixParam){
      param <- append(param, 1, after = p - 1)
    }
  }
  
  param <- as.numeric(param)
  muX <- param[1]
  muY <- param[2]
  alphaX <- param[3]
  alphaY <- param[4]
  rhoNor <- param[5]
  rhoAbn2 <- param[6]
  lengthZetaX <- dim(FPCounts)[1] - 1
  zetaX <- param[7:(6+lengthZetaX)]
  zetaY <- param[(7+lengthZetaX):length(param)]
  # print(zetaX)
  # print(zetaY)
  rhoAbn1 <- rhoNor
  rhoAbn2 <- rhoAbn2
  rhoAbn3 <- mean(c(rhoAbn1, rhoAbn2))
  rhoAbn4 <- rhoAbn3
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaAbn1 <- rbind(c(1, rhoAbn1), c(rhoAbn1, 1))
  sigmaAbn2 <- rbind(c(1, rhoAbn2), c(rhoAbn2, 1))
  sigmaAbn3 <- rbind(c(1, rhoAbn3), c(rhoAbn3, 1))
  sigmaAbn4 <- rbind(c(1, rhoAbn4), c(rhoAbn4, 1))
  
  LLNor <- 0
  LLAbn <- 0
  nBinsX <- dim(FPCounts)[1]
  nBinsY <- dim(FPCounts)[2]
  for (bX in 1:nBinsX){
    for (bY in 1:nBinsY){
      LLNor <- LLNor + FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
      LLAbn <- LLAbn + TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
    }
  }
  LL <- LLNor + LLAbn
  
  return(LL)
}

CorCBMLLNoTrns2 <- function(param, FPCounts, TPCounts){
  muX <- param[1]
  muY <- param[2]
  alphaX <- param[3]
  alphaY <- 1
  rhoNor <- param[4]
  rhoAbn2 <- param[5]
  lengthZetaX <- dim(FPCounts)[1] - 1
  zetaX <- param[6:(5+lengthZetaX)]
  zetaY <- param[(6+lengthZetaX):length(param)]
  # print(zetaX)
  # print(zetaY)
  rhoAbn1 <- rhoNor
  rhoAbn2 <- rhoAbn2
  rhoAbn3 <- mean(c(rhoAbn1, rhoAbn2))
  rhoAbn4 <- rhoAbn3
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaAbn1 <- rbind(c(1, rhoAbn1), c(rhoAbn1, 1))
  sigmaAbn2 <- rbind(c(1, rhoAbn2), c(rhoAbn2, 1))
  sigmaAbn3 <- rbind(c(1, rhoAbn3), c(rhoAbn3, 1))
  sigmaAbn4 <- rbind(c(1, rhoAbn4), c(rhoAbn4, 1))
  
  LLNor <- 0
  LLAbn <- 0
  nBinsX <- dim(FPCounts)[1]
  nBinsY <- dim(FPCounts)[2]
  for (bX in 1:nBinsX){
    for (bY in 1:nBinsY){
      LLNor <- LLNor + FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
      LLAbn <- LLAbn + TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
    }
  }
  LL <- LLNor + LLAbn
  
  return(LL)
}


CorCBMNLLGradTest <- function(muXFwd, muYFwd, alphaXFwd, alphaYFwd, rhoNorFwd, 
                              rhoAbn2Fwd, #FPCounts, TPCounts,
                              zetaXFwd1, zetaXFwd2, zetaXFwd3, zetaXFwd4, 
                              zetaYFwd1, zetaYFwd2, zetaYFwd3, zetaYFwd4
){
  allParameters <- names(formals())
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaXFwd <- unlist(mget(allParameters[which(zetaXPos == 1)]))
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaYFwd <- unlist(mget(allParameters[which(zetaYPos == 1)]))
  
  muX <- InverseValue(muXFwd, muMin, muMax)
  muY <- InverseValue(muYFwd, muMin, muMax)
  alphaX <- InverseValue(alphaXFwd, alphaMin, alphaMax)
  alphaY <- InverseValue(alphaYFwd, alphaMin, alphaMax)
  rhoNor <- InverseValue(rhoNorFwd, rhoMin, rhoMax)
  rhoAbn1 <- InverseValue(rhoNorFwd, rhoMin, rhoMax)
  rhoAbn2 <- InverseValue(rhoAbn2Fwd, rhoMin, rhoMax)
  rhoAbn3 <- mean(c(rhoAbn1, rhoAbn2))
  rhoAbn4 <- rhoAbn3
  zetaX <- InverseZetas(zetaXFwd)
  zetaY <- InverseZetas(zetaYFwd)
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaAbn1 <- rbind(c(1, rhoAbn1), c(rhoAbn1, 1))
  sigmaAbn2 <- rbind(c(1, rhoAbn2), c(rhoAbn2, 1))
  sigmaAbn3 <- rbind(c(1, rhoAbn3), c(rhoAbn3, 1))
  sigmaAbn4 <- rbind(c(1, rhoAbn4), c(rhoAbn4, 1))
  
  grad <- rep(0, length(allParameters))
  nBinsX <- dim(FPCounts)[1]
  nBinsY <- dim(FPCounts)[2]
  for (bX in 1:nBinsX){
    for (bY in 1:nBinsY){
      # LLGradNor <- LLNor + FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
      # LLGradAbn <- LLAbn + TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
      #                                               (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
      #                                               (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
      #                                               alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
      grad[1] <- grad[1] + TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                 (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                 (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                 alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (alphaX) * (1 - alphaY) * ( 
            (-dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY + 1] - rhoAbn3 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn3^2))) - 
              (-dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY + 1] - rhoAbn3 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn3^2))) - 
              (-dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY] - rhoAbn3 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn3^2))) +
              (-dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY] - rhoAbn3 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn3^2)))
          ) + 
            alphaX * alphaY * (
              (-dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY + 1] - muY - rhoAbn2 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn2^2))) - 
                (-dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY + 1] - muY - rhoAbn2 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn2^2))) - 
                (-dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY] - muY - rhoAbn2 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn2^2))) +
                (-dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY] - muY - rhoAbn2 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn2^2)))
            )
        )
      
      grad[2] <- grad[2] + TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                 (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                 (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                 alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (1 - alphaX) * alphaY * ( 
            (-dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX + 1] - rhoAbn4 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn4^2))) - 
              (-dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX + 1] - rhoAbn4 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn4^2))) - 
              (-dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX] - rhoAbn4 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn4^2))) +
              (-dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX] - rhoAbn4 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn4^2)))
          ) + 
            alphaX * alphaY * (
              (-dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX + 1] - muX - rhoAbn2 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn2^2))) - 
                (-dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX + 1] - muX - rhoAbn2 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn2^2))) - 
                (-dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX] - muX - rhoAbn2 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn2^2))) +
                (-dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX] - muX - rhoAbn2 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn2^2)))
            )
        )
      
      grad[3] <- grad[3] + TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                 (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                 (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                 alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (-(1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
             (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
             - (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
             alphaY * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
        )
      
      grad[4] <- grad[4] + TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                 (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                 (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                 alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (-(1 - alphaX) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
             (1 - alphaX) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
             - (alphaX) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
             alphaX * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
        )
      
      grad[5] <- grad[5] + 
        FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
        (
          dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor) - 
            dmvnorm(c(zetaX[bX + 1], zetaY[bY]), sigma = sigmaNor) - 
            dmvnorm(c(zetaX[bX], zetaY[bY + 1]), sigma = sigmaNor) + 
            dmvnorm(c(zetaX[bX], zetaY[bY]), sigma = sigmaNor)
          
        ) + 
        TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                              (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                              (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                              alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (1 - alphaX) * (1 - alphaY) * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) -
                                           dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) -
                                           dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(0, 0), sigma = sigmaAbn1) +
                                           dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(0, 0), sigma = sigmaAbn1)) + 
            (alphaX) * (1 - alphaY) * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) - 
                                         dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) -
                                         dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                         dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3)) / 2+ 
            (1 - alphaX) * (alphaY) * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) - 
                                         dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) -
                                         dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4) +
                                         dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4)) / 2
        )
      
      grad[6] <- grad[6] + 
        TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                              (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                              (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                              alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (alphaX) * (1 - alphaY) * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) - 
                                       dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) -
                                       dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                       dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3)) / 2 + 
            (1 - alphaX) * (alphaY) * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) - 
                                         dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) -
                                         dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4) +
                                         dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4)) / 2 + 
            alphaX * alphaY * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2) -
                                 dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2) -
                                 dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(muX, muY), sigma = sigmaAbn2) +
                                 dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(muX, muY), sigma = sigmaAbn2))
        )
      
      if (bX == 1){
        grad[7] <- grad[7] + 
          FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
          (
            dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor) - 
              dmvnorm(c(zetaX[bX + 1], zetaY[bY]), sigma = sigmaNor) - 
              dmvnorm(c(zetaX[bX], zetaY[bY + 1]), sigma = sigmaNor) + 
              dmvnorm(c(zetaX[bX], zetaY[bY]), sigma = sigmaNor)
            
          ) + 
          TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
          (
            (1 - alphaX) * (1 - alphaY) * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) -
                                             dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) -
                                             dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(0, 0), sigma = sigmaAbn1) +
                                             dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(0, 0), sigma = sigmaAbn1)) + 
              (alphaX) * (1 - alphaY) * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) - 
                                           dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) -
                                           dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                           dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3)) / 2+ 
              (1 - alphaX) * (alphaY) * (dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) - 
                                           dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) -
                                           dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4) +
                                           dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4)) / 2
          )
      }
    }
  }
  
  return(-grad)
}

# muX <- muIniX
# muY <- muIniY
# alpha <- alphaIni
# zetaX <- zetaIniX
# zetaY <- zetaIniY
# rhoNor <- rhoIniNor
# rhoAbn <- rhoIniAbn
# sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
# sigmaAbn <- rbind(c(1, rhoAbn), c(rhoAbn, 1))

# CBMAbnCdf <- function(muX, muY, alpha, sigmaNor, sigmaAbn, zetaXLower, zetaXUpper, zetaYLower, zetaYUpper){
#   int <- adaptIntegrate(CBMAbnPdf, c(zetaXLower, zetaYLower), c(zetaXUpper, zetaYUpper), muX = muX, muY = muY, alpha = alpha, sigmaNor = sigmaNor, sigmaAbn = sigmaAbn, tol = 1e-6)
#   return(int$integral)
# }
# 
# CBMAbnIntInfUpper <- function(txy, zetaXLower, zetaYLower, muX, muY, alpha, sigmaNor, sigmaAbn){
#   tx <- txy[1]
#   ty <- txy[2]
#   x <- tx / (1 - tx) + zetaXLower
#   y <- ty / (1 - ty) + zetaYLower
#   xy <- c(x, y)
#   abnPdf <- CBMAbnPdf(xy, muX, muY, alpha, sigmaNor, sigmaAbn) / (1 - tx)^2 / (1 - ty)^2
#   return(abnPdf)
# }
# 
# CBMAbnIntInfLower <- function(txy, zetaXUpper, zetaYUpper, muX, muY, alpha, sigmaNor, sigmaAbn){
#   tx <- txy[1]
#   ty <- txy[2]
#   x <- tx / (1 + tx) + zetaXUpper
#   y <- ty / (1 + ty) + zetaYUpper
#   xy <- c(x, y)
#   abnPdf <- CBMAbnPdf(xy, muX, muY, alpha, sigmaNor, sigmaAbn) / (1 + tx)^2 / (1 + ty)^2
#   return(abnPdf)
# }
# 
# CBMAbnIntInf <- function(txy, muX, muY, alpha, sigmaNor, sigmaAbn){
#   tx <- txy[1]
#   ty <- txy[2]
#   x <- tx / (1 - tx^2)
#   y <- ty / (1 - ty^2)
#   xy <- c(x, y)
#   abnPdf <- CBMAbnPdf(xy, muX, muY, alpha, sigmaNor, sigmaAbn) * (1 + tx^2) / (1 - tx^2)^2 * (1 + ty^2) / (1 - ty^2)^2
#   return(abnPdf)
# }
# 
# CBMAbnPdf <- function(xy){
#   x <- xy[1]; y <- xy[2]
#   abnPdf <- dmvnorm(c(x, y))
#   return(abnPdf)
# }




# CBMOuterInt <- function(muX, muY, alphaX, alphaY, rhoAbn, zeta1X, zeta2X, zeta1Y, zeta2Y){
#   intgrl <- integrate(CBMInnerFun, lower = zeta1X, upper = zeta2X, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoAbn = rhoAbn, zeta1Y = zeta1Y, zeta2Y = zeta2Y, rel.tol = .Machine$double.eps^0.5)$value
#   # intgrl <- integral(CBMInnerFun, xmin = zeta1X, xmax = zeta2X, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoAbn = rhoAbn, zeta1Y = zeta1Y, zeta2Y = zeta2Y)
#   return(intgrl)
# }
# 
# CBMInnerFun <- function(x, muX, muY, alphaX, alphaY, rhoAbn, zeta1Y, zeta2Y){
#   return(sapply(x, CBMInnerInt, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoAbn = rhoAbn, zeta1Y = zeta1Y, zeta2Y = zeta2Y))
# }
# 
# CBMInnerInt <- function(x, muX, muY, alphaX, alphaY, rhoAbn, zeta1Y, zeta2Y){
#   intgrl <- integrate(CBMAbnFun, lower = zeta1Y, upper = zeta2Y, x = x, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoAbn = rhoAbn, rel.tol = .Machine$double.eps^0.5)$value
#   # intgrl <- integral(CBMAbnFun, xmin = zeta1Y, xmax = zeta2Y, x = x, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoAbn = rhoAbn)
#   return(intgrl)
# }
# 
# CBMAbnFun <- function(x, y, muX, muY, alphaX, alphaY, rhoAbn){
#   return(sapply(y, CBMAbnPdf, x = x, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoAbn = rhoAbn))
# }
# 
# CBMCondPdf <- function(x, y, alphaX, muX, muY, rhoNor, rhoMu){
#   muAbnX <- (1 - alphaX) * rhoNor * y + alphaX * (muX + rhoMu * (y - muY))
#   sigmaAbnX <- sqrt((1 - alphaX)^2 * (1 - rhoNor^2) + alphaX^2 * (1 - rhoMu^2))
#   condPdf <- dnorm(x, muAbnX, sigmaAbnX)
#   return(condPdf)
# }