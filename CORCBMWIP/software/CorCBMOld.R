CorCBM <- function(FPCounts, TPCounts, aucX, aucY, K1, K2, rhoIniNor, rhoIniAbn){
  iniParamX <- IniParam(FPCounts[1, ], TPCounts[1, ], aucX, K1, K2)
  alphaIniX <- iniParamX$alphaIni
  muIniX <- iniParamX$muIni
  zetaIniX <- iniParamX$zetaIni
  
  iniParamY <- IniParam(FPCounts[2, ], FPCounts[2, ], aucY, K1, K2)
  alphaIniY <- iniParamY$alphaIni
  muIniY <- iniParamY$muIni
  zetaIniY <- iniParamY$zetaIni
  
  rhoIniMu <- rhoIniAbn
}

IniParam <- function(fpCounts, tpCounts, auc, K1, K2){
  alphaIni <- 1
  muIni <- qnorm(auc) * sqrt(2)
  fpCountsCum <- cumsum(fpCounts)[1:(length(fpCounts) - 1)]
  tpCountsCum <- cumsum(tpCounts)[1:(length(tpCounts) - 1)]
  zetaNorIni <- qnorm(fpCountsCum / K1)
  zetaAbnIni <- qnorm(tpCountsCum / K2, mean = muIni)
  for (z in 1:length(zetaNorIni)){
    if (z == 1){
      zetaNorIni[z] <- max(zetaNorIni[z], -20)
      zetaAbnIni[z] <- max(zetaAbnIni[z], -20)
    }else{
      zetaNorIni[z] <- max(zetaNorIni[z], zetaNorIni[z - 1] + 0.001)
      zetaAbnIni[z] <- max(zetaAbnIni[z], zetaAbnIni[z - 1] + 0.001)
    }
  }
  zetaIni <- (zetaNorIni + zetaAbnIni) / 2
  return(list(alphaIni = alphaIni,
              muIni = muIni,
              zetaIni = zetaIni))
}

muX <- muIniX
muY <- muIniY
alphaX <- alphaIniX
alphaY <- alphaIniY
zetaX <- zetaIniX
zetaY <- zetaIniY
rhoNor <- rhoIniNor
rhoAbn <- rhoIniAbn
rhoMu <- rhoIniMu

CorCBMNLL <- function(muX, muY, alphaX, alphaY, rhoNor, rhoAbn, zetaX, zetaY, FPCounts, TPCounts){
  nBins <- length(FPCounts[1, ])
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaMu <- rbind(c(1, rhoAbn), c(rhoAbn, 1))
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  LL1 <- 0
  LL2 <- 0
  bX <- 1
  bY <- 1
  for (bX in 1:nBins){
    for (bY in 1:nBins){
      LLNor <- (FPCounts[1, bX] + FPCounts[2, bY]) * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
      LL2Int <- 
      LL2
    }
  }
}

CBMOuterInt <- function(muX, muY, alphaX, alphaY, rhoNor, rhoAbn, rhoMu, zeta1X, zeta2X, zeta1Y, zeta2Y){
  # intgrl <- integrate(CBMInnerFun, lower = zeta1X, upper = zeta2X, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoNor = rhoNor, rhoAbn = rhoAbn, rhoMu = rhoMu, zeta1Y = zeta1Y, zeta2Y = zeta2Y, rel.tol = .Machine$double.eps^0.5)$value
  intgrl <- integral(CBMInnerFun, xmin = zeta1X, xmax = zeta2X, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoNor = rhoNor, rhoAbn = rhoAbn, rhoMu = rhoMu, zeta1Y = zeta1Y, zeta2Y = zeta2Y)
  return(intgrl)
}

CBMInnerFun <- function(x, muX, muY, alphaX, alphaY, rhoNor, rhoAbn, rhoMu, zeta1Y, zeta2Y){
  return(sapply(x, CBMInnerInt, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoNor = rhoNor, rhoAbn = rhoAbn, rhoMu = rhoMu, zeta1Y = zeta1Y, zeta2Y = zeta2Y))
}

CBMInnerInt <- function(x, muX, muY, alphaX, alphaY, rhoNor, rhoAbn, rhoMu, zeta1Y, zeta2Y){
  # intgrl <- integrate(CBMAbnFun, lower = zeta1Y, upper = zeta2Y, x = x, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoNor = rhoNor, rhoAbn = rhoAbn, rhoMu = rhoMu, rel.tol = .Machine$double.eps^0.5)$value
  intgrl <- integral(CBMAbnFun, xmin = zeta1Y, xmax = zeta2Y, x = x, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoNor = rhoNor, rhoAbn = rhoAbn, rhoMu = rhoMu)
  return(intgrl)
}

CBMOuterFun <- function(x, y, muX, muY, alphaX, alphaY, rhoNor, rhoAbn, rhoMu){
  return(mapply(CBMAbnPdf, x, y, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoNor = rhoNor, rhoAbn = rhoAbn, rhoMu = rhoMu))
}

CBMAbnFun <- function(x, y, muX, muY, alphaX, alphaY, rhoNor, rhoAbn, rhoMu){
  return(sapply(y, CBMAbnPdf, x = x, muX = muX, muY = muY, alphaX = alphaX, alphaY = alphaY, rhoNor = rhoNor, rhoAbn = rhoAbn, rhoMu = rhoMu))
}

CBMAbnPdf <- function(x, y, muX, muY, alphaX, alphaY, rhoNor, rhoAbn, rhoMu){
  muAbnX <- (1 - alphaX) * rhoNor * y + alphaX * (muX + rhoMu * (y - muY))
  sigmaAbnX <- sqrt((1 - alphaX)^2 * (1 - rhoNor^2) + alphaX^2 * (1 - rhoMu^2))
  
  muAbnY <- (1 - alphaY) * rhoNor * x + alphaY * (muY + rhoMu * (x - muX))
  sigmaAbnY <- sqrt((1 - alphaY)^2 * (1 - rhoNor^2) + alphaY^2 * (1 - rhoMu^2))
  
  sigmaAbn <- rbind(c(sigmaAbnX^2, rhoAbn * sigmaAbnX * sigmaAbnY), c(rhoAbn * sigmaAbnX * sigmaAbnY, sigmaAbnY^2))
  
  abnPdf <- dmvnorm(c(x, y), c(muAbnX, muAbnY), sigmaAbn)
  return(abnPdf)
}

testFun <- function(x, y){
  return(mapply(testFun2, x, y))
}

testFun2 <- function(x,y){
  return(dmvnorm(c(x, y), mean = c(1, 2)))
}

CBMCondPdf <- function(x, y, alphaX, muX, muY, rhoNor, rhoMu){
  muAbnX <- (1 - alphaX) * rhoNor * y + alphaX * (muX + rhoMu * (y - muY))
  sigmaAbnX <- sqrt((1 - alphaX)^2 * (1 - rhoNor^2) + alphaX^2 * (1 - rhoMu^2))
  condPdf <- dnorm(x, muAbnX, sigmaAbnX)
  return(condPdf)
}