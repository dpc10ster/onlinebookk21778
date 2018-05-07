CorCBMNR <- function(FPCounts, TPCounts, aucX, aucY, K1, K2, rhoIniNor, rhoIniAbn){
  iniParamX <- IniParam(FPCounts[1, ], TPCounts[1, ], aucX, K1, K2)
  muIniX <- iniParamX$muIni
  zetaIniX <- iniParamX$zetaIni
  
  iniParamY <- IniParam(FPCounts[2, ], TPCounts[2, ], aucY, K1, K2)
  muIniY <- iniParamY$muIni
  zetaIniY <- iniParamY$zetaIni
  alphaIni <- 0.9
  
  parameters <- c(list(muIniX, muIniY, alphaIni, rhoIniNor, rhoIniAbn), as.list(zetaIniX), as.list(zetaIniY))
  namesVector <- c("muX", "muY", "alpha", "rhoNor", "rhoAbn")
  for (z in 1:length(zetaIniX)){
    namesVector <- c(namesVector, paste0("zetaX", z))
  }
  for (z in 1:length(zetaIniY)){
    namesVector <- c(namesVector, paste0("zetaY", z))
  }
  names(parameters) <- namesVector
  
  CorCBMNLLNew <- AddZetaX(CorCBMNLL, length(zetaIniX))
  CorCBMNLLNew <- AddZetaY(CorCBMNLLNew, length(zetaIniY))
  arguments <- c(parameters, list(FPCounts = FPCounts), list(TPCounts = TPCounts))
  do.call(CorCBMNLLNew, arguments)
  nBins <- dim(FPCounts)[2]
  FPPrEst <- rep(NA, nBins^2)
  FPComb <- FPPrEst
  TPPrEst <- FPPrEst
  TPComb <- FPComb
  bInx <- 1
  for (bX in 1:nBins){
    for (bY in 1:nBins){
      FPComb[bInx] <- FPCounts[1, bX] + FPCounts[2, bY]
      TPComb[bInx] <- TPCounts[1, bX] + TPCounts[2, bY]
      bInx <- bInx + 1
    }
  }
  Kdd1 <- sum(FPComb)
  Kdd2 <- sum(TPComb)
  
  allCoef <- ret@coef
  muX <- as.numeric(InverseValue(allCoef[1], muMin, muMax))
  muY <- as.numeric(InverseValue(allCoef[2], muMin, muMax))
  alpha <- as.numeric(InverseValue(allCoef[3], alphaMin, alphaMax))
  rhoNor <- as.numeric(InverseValue(allCoef[4], rhoMin, rhoMax))
  rhoAbn <- as.numeric(InverseValue(allCoef[5], rhoMin, rhoMax))
  zetaX <- as.numeric(InverseZetas(allCoef[6:(5 + length(zetaIniX))]))
  zetaY <- as.numeric(InverseZetas(allCoef[(6 + length(zetaIniX)):length(allCoef)]))
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaAbn <- rbind(c(1, rhoAbn), c(rhoAbn, 1))
  
  
  plot(FPPrEmp, FPPrEst)
}

IniParam <- function(fpCounts, tpCounts, auc, K1, K2){
  muIni <- qnorm(auc) * sqrt(2)
  fpCountsCum <- cumsum(fpCounts)[1:(length(fpCounts) - 1)]
  fpCountsCum[fpCountsCum == 0] <- pnorm(-zetaMin)
  tpCountsCum <- cumsum(tpCounts)[1:(length(tpCounts) - 1)]
  tpCountsCum[tpCountsCum == 0] <- pnorm(-zetaMin, mean = muIni)
  
  zetaNorIni <- qnorm(fpCountsCum / K1)
  zetaAbnIni <- qnorm(tpCountsCum / K2, mean = muIni)
  for (z in 1:length(zetaNorIni)){
    zetaNorIni[z] <- min(max(zetaNorIni[z], zetaNorIni[z - 1] + 0.001), zetaMax)
    zetaAbnIni[z] <- min(max(zetaAbnIni[z], zetaAbnIni[z - 1] + 0.001), zetaMax)
  }
  zetaIni <- (zetaNorIni + zetaAbnIni) / 2
  # zetaIni <- zetaNorIni
  return(list(muIni = muIni,
              zetaIni = zetaIni))
}


CorCBMNLLTest <- function(muXF, muY, alpha, rhoNor, rhoAbn, FPCounts, TPCounts,
                          zetaX1, zetaX2, zetaX3, zetaX4, zetaY1, 
                          zetaY2, zetaY3, zetaY4){
  allParameters <- names(formals())
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(mget(allParameters[which(zetaXPos == 1)]))
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(mget(allParameters[which(zetaYPos == 1)]))
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaAbn <- rbind(c(1, rhoAbn), c(rhoAbn, 1))
  
  LLNor <- 0
  LLAbn <- 0
  nBins <- dim(FPCounts)[2]
  for (bX in 1:nBins){
    for (bY in 1:nBins){
      LLNor <- LLNor + (FPCounts[1, bX] + FPCounts[2, bY]) * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
      LLAbn <- LLAbn + (TPCounts[1, bX] + TPCounts[2, bY]) * log((1 - alpha) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaAbn) + 
                                                                   alpha * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn))
    }
  }
  LL <- LLNor + LLAbn
  if (is.infinite(LL)){
    dummyStop <- 1
  }
  return(-LL)
}

CorCBMNLL <- function(muX, muY, alpha, rhoNor, rhoAbn, FPCounts, TPCounts){
  allParameters <- names(formals())
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(mget(allParameters[which(zetaXPos == 1)]))
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(mget(allParameters[which(zetaYPos == 1)]))
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaAbn <- rbind(c(1, rhoAbn), c(rhoAbn, 1))
  
  LLNor <- 0
  LLAbn <- 0
  nBins <- dim(FPCounts)[2]
  for (bX in 1:nBins){
    for (bY in 1:nBins){
      LLNor <- LLNor + (FPCounts[1, bX] + FPCounts[2, bY]) * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
      LLAbn <- LLAbn + (TPCounts[1, bX] + TPCounts[2, bY]) * log((1 - alpha) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaAbn) + 
                                                                   alpha * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn))
    }
  }
  LL <- LLNor + LLAbn
  return(-LL)
}

