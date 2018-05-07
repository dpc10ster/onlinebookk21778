CorCBM <- function(FPCounts, TPCounts, aucX, aucY, K1, K2, rhoIniNor, rhoIniAbn2){
  iniParamX <- IniParam(rowSums(FPCounts), rowSums(TPCounts), aucX, K1, K2)
  muIniX <- iniParamX$muIni
  zetaIniX <- iniParamX$zetaIni
  alphaIniX <- iniParamX$alphaIni
  # zetaIniX <- c(zetaIniX[1], zetaIniX[2:4] - zetaIniX[1:3])
  
  iniParamY <- IniParam(colSums(FPCounts), colSums(TPCounts), aucY, K1, K2)
  muIniY <- iniParamY$muIni
  zetaIniY <- iniParamY$zetaIni
  alphaIniY <- iniParamY$alphaIni
  # zetaIniY <- c(zetaIniY[1], zetaIniY[2:4] - zetaIniY[1:3])
  
  parameters <- c(list(muIniX, muIniY, alphaIniX, alphaIniY, rhoIniNor, rhoIniAbn2), as.list(zetaIniX), as.list(zetaIniY))
  namesVector <- c("muX", "muY", "alphaX", "alphaY","rhoNor", "rhoAbn2")
  for (z in 1:length(zetaIniX)){
    namesVector <- c(namesVector, paste0("zetaX", z))
  }
  for (z in 1:length(zetaIniY)){
    namesVector <- c(namesVector, paste0("zetaY", z))
  }
  names(parameters) <- namesVector
  
  # CorCBMNLLNew <- AddZetaXFwd(CorCBMNLL, length(zetaIniX))
  # CorCBMNLLNew <- AddZetaYFwd(CorCBMNLLNew, length(zetaIniY))
  ret <- mle2(CorCBMNLLTest, start = parameters, method = "BFGS", data = list(FPCounts = FPCounts, TPCounts = TPCounts), gr = CorCBMNLLGradTest,
              control = list(parscale = c(rep(1, 6), rep(0.01, 8))))
  # deltaZeta <- 1e-3
  # ret <- mle2(CorCBMNLLTest, start = parameters, method = "L-BFGS-B", data = list(FPCounts = FPCounts, TPCounts = TPCounts), #gr = CorCBMNLLGradTest, 
  #             lower = c(muMin, muMin, alphaMin, alphaMin, rhoMin, rhoMin, zetaMin, deltaZeta, deltaZeta, deltaZeta, zetaMin, deltaZeta, deltaZeta, deltaZeta), 
  #             upper = c(muMax, muMax, alphaMax, alphaMax, rhoMax, rhoMax, zetaMax, zetaMax, zetaMax, zetaMax, zetaMax, zetaMax, zetaMax, zetaMax))
  # ret <- mle2(CorCBMNLLNew, start = parameters, method = "BFGS", data = list(FPCounts = FPCounts, TPCounts = TPCounts))
  
  allCoef <- ret@coef
  muX <- as.numeric((allCoef[1]))
  muY <- as.numeric((allCoef[2]))
  alphaX <- as.numeric((allCoef[3]))
  alphaY <- as.numeric((allCoef[4]))
  rhoNor <- as.numeric((allCoef[5]))
  rhoAbn2 <- as.numeric((allCoef[6]))
  zetaX <- as.numeric((allCoef[7:(6 + length(zetaIniX))]))
  zetaY <- as.numeric((allCoef[(7 + length(zetaIniX)):length(allCoef)]))
  print(c(muX, muY, alphaX, alphaY, rhoNor, rhoAbn2, zetaX, zetaY))
  
  
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  rhoAbn1 <- rhoNor
  rhoAbn3 <- mean(c(rhoAbn1, rhoAbn2))
  rhoAbn4 <- rhoAbn3
  sigmaAbn1 <- rbind(c(1, rhoAbn1), c(rhoAbn1, 1))
  sigmaAbn2 <- rbind(c(1, rhoAbn2), c(rhoAbn2, 1))
  sigmaAbn3 <- rbind(c(1, rhoAbn3), c(rhoAbn3, 1))
  sigmaAbn4 <- rbind(c(1, rhoAbn4), c(rhoAbn4, 1))
  
  nBins <- dim(FPCounts)[2]
  FPPrEst <- rep(NA, nBins^2)
  FPNum <- FPPrEst
  TPPrEst <- FPPrEst
  TPNum <- FPNum
  bInx <- 1
  for (bX in 1:nBins){
    for (bY in 1:nBins){
      FPPrEst[bInx] <- pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)
      FPNum[bInx] <- FPCounts[bX, bY]
      TPPrEst[bInx] <- (1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
        (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
        (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
        alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)
      TPNum[bInx] <- TPCounts[bX, bY]
      bInx <- bInx + 1
    }
  }
  FPPrEmp <- FPNum / sum(FPNum)
  TPPrEmp <- TPNum / sum(TPNum)
  plot(FPPrEmp, FPPrEst)
  lines(c(0, max(c(FPPrEmp, FPPrEst))), c(0, max(c(FPPrEmp, FPPrEst))))
  plot(TPPrEmp, TPPrEst)
  lines(c(0, max(c(TPPrEmp, TPPrEst))), c(0, max(c(TPPrEmp, TPPrEst))))
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
  tpf <- cumsum(rev(tpCounts))/sum(tpCounts)
  alphaIni <- mean(rev(tpf)[2:3])
  
  muIniFwd <- ForwardValue(muIni, muMin, muMax)
  alphaIniFwd <- ForwardValue(alphaIni, alphaMin, alphaMax)
  zetaIniFwd <- ForwardZetas(zetaIni)
  parameters <- c(list(muIniFwd, alphaIniFwd), as.list(zetaIniFwd))
  namesVector <- c("muFwd", "alphaFwd")
  for (z in 1:length(zetaIni)){
    namesVector <- c(namesVector, paste0("zetaFwd", z))
  }
  names(parameters) <- namesVector
  ret <- mle2(CBMNLLTest, start = parameters, method = "BFGS", data = list(FPCounts = fpCounts, TPCounts = tpCounts))
  allCoef <- ret@coef
  muFwd <- allCoef[1]
  alphaFwd <- allCoef[2]
  zetaFwd <- allCoef[3:length(allCoef)]
  
  muIni <- InverseValue(muFwd, muMin, muMax)
  alphaIni <- InverseValue(alphaFwd, alphaMin, alphaMax)
  zetaIni <- InverseZetas(zetaFwd)
  
  return(list(muIni = muIni,
              zetaIni = zetaIni, 
              alphaIni = alphaIni))
}


CorCBMNLLTest <- function(muX, muY, alphaX, alphaY, rhoNor, 
                          rhoAbn2, #FPCounts, TPCounts,
                          zetaX1, zetaX2, zetaX3, zetaX4, 
                          zetaY1, zetaY2, zetaY3, zetaY4
){
  allParameters <- names(formals())
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(mget(allParameters[which(zetaXPos == 1)]))
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(mget(allParameters[which(zetaYPos == 1)]))
  
  rhoAbn1 <- rhoNor
  rhoAbn2 <- rhoAbn2
  rhoAbn3 <- mean(c(rhoAbn1, rhoAbn2))
  rhoAbn4 <- rhoAbn3
  
  # zetaX <- cumsum(zetaX)
  # zetaY <- cumsum(zetaY)
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  sigmaNor <- rbind(c(1, rhoNor), c(rhoNor, 1))
  sigmaAbn1 <- rbind(c(1, rhoAbn1), c(rhoAbn1, 1))
  sigmaAbn2 <- rbind(c(1, rhoAbn2), c(rhoAbn2, 1))
  sigmaAbn3 <- rbind(c(1, rhoAbn3), c(rhoAbn3, 1))
  sigmaAbn4 <- rbind(c(1, rhoAbn4), c(rhoAbn4, 1))
  
  LLNor <- 0
  LLAbn <- 0
  nBins <- dim(FPCounts)[2]
  for (bX in 1:nBins){
    for (bY in 1:nBins){
      if (bX == 5 && bY == 5){
        breakPoint <- 1
      }
      LLNor <- LLNor + FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
      LLAbn <- LLAbn + TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
    }
  }
  LL <- LLNor + LLAbn
  
  # if (is.infinite(LL)){
  #   LLNor <- 0
  #   LLAbn <- 0
  #   for (bX in 1:nBins){
  #     for (bY in 1:nBins){
  #       LLNor <- LLNor + FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))
  #       LLAbn <- LLAbn + TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
  #                                                 (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
  #                                                 (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
  #                                                 alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
  #       if (any(is.infinite(c(LLNor, LLAbn)))){
  #         dummyStop <- 1
  #       }
  #     }
  #   }
  # }
  return(-LL)
}


CBMNLLTest <- function(muFwd, alphaFwd, FPCounts, TPCounts,
                       zetaFwd1, zetaFwd2, zetaFwd3, zetaFwd4
){
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



CorCBMNLLGradTest <- function(muX, muY, alphaX, alphaY, rhoNor, 
                              rhoAbn2, #FPCounts, TPCounts,
                              zetaX1, zetaX2, zetaX3, zetaX4, 
                              zetaY1, zetaY2, zetaY3, zetaY4
){
  allParameters <- names(formals())
  zetaXPos <- regexpr("zetaX", allParameters)
  zetaX <- unlist(mget(allParameters[which(zetaXPos == 1)]))
  zetaYPos <- regexpr("zetaY", allParameters)
  zetaY <- unlist(mget(allParameters[which(zetaYPos == 1)]))
  
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
  
  grad <- rep(0, length(allParameters))
  nBins <- dim(FPCounts)[2]
  for (bX in 1:nBins){
    for (bY in 1:nBins){
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
            ifelse(zetaX[bX + 1] %in% c(-Inf, Inf), 0, (-dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY + 1] - rhoAbn3 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn3^2)))) - 
              ifelse(zetaX[bX] %in% c(-Inf, Inf), 0, (-dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY + 1] - rhoAbn3 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn3^2)))) - 
              ifelse(zetaX[bX + 1] %in% c(-Inf, Inf), 0, (-dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY] - rhoAbn3 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn3^2)))) +
              ifelse(zetaX[bX] %in% c(-Inf, Inf), 0, (-dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY] - rhoAbn3 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn3^2))))
          ) + 
            alphaX * alphaY * (
              ifelse(zetaX[bX + 1] %in% c(-Inf, Inf), 0, (-dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY + 1] - muY - rhoAbn2 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn2^2)))) - 
                ifelse(zetaX[bX] %in% c(-Inf, Inf), 0, (-dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY + 1] - muY - rhoAbn2 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn2^2)))) - 
                ifelse(zetaX[bX + 1] %in% c(-Inf, Inf), 0, (-dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY] - muY - rhoAbn2 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn2^2)))) +
                ifelse(zetaX[bX] %in% c(-Inf, Inf), 0, (-dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY] - muY - rhoAbn2 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn2^2))))
            )
        )
      
      grad[2] <- grad[2] + TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                                 (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                                 (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                                 alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (1 - alphaX) * alphaY * ( 
            ifelse(zetaY[bY + 1] %in% c(-Inf, Inf), 0, (-dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX + 1] - rhoAbn4 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn4^2)))) - 
              ifelse(zetaY[bY] %in% c(-Inf, Inf), 0, (-dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX + 1] - rhoAbn4 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn4^2)))) - 
              ifelse(zetaY[bY + 1] %in% c(-Inf, Inf), 0, (-dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX] - rhoAbn4 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn4^2)))) +
              ifelse(zetaY[bY] %in% c(-Inf, Inf), 0, (-dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX] - rhoAbn4 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn4^2))))
          ) + 
            alphaX * alphaY * (
              ifelse(zetaY[bY + 1] %in% c(-Inf, Inf), 0, (-dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX + 1] - muX - rhoAbn2 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn2^2)))) - 
                ifelse(zetaY[bY] %in% c(-Inf, Inf), 0, (-dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX + 1] - muX - rhoAbn2 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn2^2)))) - 
                ifelse(zetaY[bY + 1] %in% c(-Inf, Inf), 0, (-dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX] - muX - rhoAbn2 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn2^2)))) +
                ifelse(zetaY[bY] %in% c(-Inf, Inf), 0, (-dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX] - muX - rhoAbn2 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn2^2))))
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
             - (alphaX) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
             (1 - alphaX) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
             alphaX * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
        )
      
      grad[5] <- grad[5] + 
        FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
        (
          Dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor) - 
            Dmvnorm(c(zetaX[bX + 1], zetaY[bY]), sigma = sigmaNor) - 
            Dmvnorm(c(zetaX[bX], zetaY[bY + 1]), sigma = sigmaNor) + 
            Dmvnorm(c(zetaX[bX], zetaY[bY]), sigma = sigmaNor)
          
        ) + 
        TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                              (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                              (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                              alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (1 - alphaX) * (1 - alphaY) * (Dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) -
                                           Dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) -
                                           Dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(0, 0), sigma = sigmaAbn1) +
                                           Dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(0, 0), sigma = sigmaAbn1)) + 
            (alphaX) * (1 - alphaY) * (Dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) - 
                                         Dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) -
                                         Dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                         Dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3)) / 2+ 
            (1 - alphaX) * (alphaY) * (Dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) - 
                                         Dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) -
                                         Dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4) +
                                         Dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4)) / 2
        )
      
      grad[6] <- grad[6] + 
        TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                              (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                              (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                              alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        (
          (alphaX) * (1 - alphaY) * (Dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) - 
                                       Dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) -
                                       Dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                       Dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(muX, 0), sigma = sigmaAbn3)) / 2 + 
            (1 - alphaX) * (alphaY) * (Dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) - 
                                         Dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) -
                                         Dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4) +
                                         Dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(0, muY), sigma = sigmaAbn4)) / 2 + 
            alphaX * alphaY * (Dmvnorm(c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2) -
                                 Dmvnorm(c(zetaX[bX], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2) -
                                 Dmvnorm(c(zetaX[bX + 1], zetaY[bY]), mean = c(muX, muY), sigma = sigmaAbn2) +
                                 Dmvnorm(c(zetaX[bX], zetaY[bY]), mean = c(muX, muY), sigma = sigmaAbn2))
        )
      
      delta <- 1e-10
      if (bX > 1){
        grad[6 + bX - 1] <- grad[6 + bX - 1] + 
          FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
          (
            - dnorm(zetaX[bX]) * pnorm((zetaY[bY + 1] - rhoNor * zetaX[bX]) / sqrt(1 - rhoNor^2)) + 
              dnorm(zetaX[bX]) * pnorm((zetaY[bY] - rhoNor * zetaX[bX]) / sqrt(1 - rhoNor^2))
          ) + 
          TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
          (
            (1 - alphaX) * (1 - alphaY) * (- dnorm(zetaX[bX ]) * pnorm((zetaY[bY + 1] - rhoAbn1 * zetaX[bX]) / sqrt(1 - rhoAbn1^2)) + 
                                             dnorm(zetaX[bX]) * pnorm((zetaY[bY] - rhoAbn1 * zetaX[bX]) / sqrt(1 - rhoAbn1^2))) + 
              (alphaX) * (1 - alphaY) * (- dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY + 1] - rhoAbn3 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn3^2)) + 
                                           dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY] - rhoAbn3 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn3^2))) + 
              (1 - alphaX) * (alphaY) * (- dnorm(zetaX[bX]) * pnorm((zetaY[bY + 1] - muY - rhoAbn4 * zetaX[bX]) / sqrt(1 - rhoAbn4^2)) + 
                                           dnorm(zetaX[bX]) * pnorm((zetaY[bY] - muY - rhoAbn4 * zetaX[bX]) / sqrt(1 - rhoAbn4^2))) + 
              (alphaX) * (alphaY) * (- dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY + 1] - muY - rhoAbn2 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn2^2)) + 
                                       dnorm(zetaX[bX] - muX) * pnorm((zetaY[bY] - muY - rhoAbn2 * (zetaX[bX] - muX)) / sqrt(1 - rhoAbn2^2)))
          )
      }
      if (bX < nBins){
        grad[6 + bX] <- grad[6 + bX] + 
          FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
          (
            dnorm(zetaX[bX + 1]) * pnorm((zetaY[bY + 1] - rhoNor * zetaX[bX + 1]) / sqrt(1 - rhoNor^2)) - 
              dnorm(zetaX[bX + 1]) * pnorm((zetaY[bY] - rhoNor * zetaX[bX + 1]) / sqrt(1 - rhoNor^2))
          ) + 
          TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
          (
            (1 - alphaX) * (1 - alphaY) * (dnorm(zetaX[bX + 1]) * pnorm((zetaY[bY + 1] - rhoAbn1 * zetaX[bX + 1]) / sqrt(1 - rhoAbn1^2)) - 
                                             dnorm(zetaX[bX + 1]) * pnorm((zetaY[bY] - rhoAbn1 * zetaX[bX + 1]) / sqrt(1 - rhoAbn1^2))) + 
              (alphaX) * (1 - alphaY) * (dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY + 1] - rhoAbn3 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn3^2)) - 
                                           dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY] - rhoAbn3 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn3^2))) + 
              (1 - alphaX) * (alphaY) * (dnorm(zetaX[bX + 1]) * pnorm((zetaY[bY + 1] - muY - rhoAbn4 * zetaX[bX + 1]) / sqrt(1 - rhoAbn4^2)) - 
                                           dnorm(zetaX[bX + 1]) * pnorm((zetaY[bY] - muY - rhoAbn4 * zetaX[bX + 1]) / sqrt(1 - rhoAbn4^2))) + 
              (alphaX) * (alphaY) * (dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY + 1] - muY - rhoAbn2 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn2^2)) - 
                                       dnorm(zetaX[bX + 1] - muX) * pnorm((zetaY[bY] - muY - rhoAbn2 * (zetaX[bX + 1] - muX)) / sqrt(1 - rhoAbn2^2)))
          )
      }
        
      if (bY > 1){
        grad[10 + bY - 1] <- grad[10 + bY - 1] + 
          FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
          (
            - dnorm(zetaY[bY]) * pnorm((zetaX[bX + 1] - rhoNor * zetaY[bY]) / sqrt(1 - rhoNor^2)) + 
              dnorm(zetaY[bY]) * pnorm((zetaX[bX] - rhoNor * zetaY[bY]) / sqrt(1 - rhoNor^2))
          ) + 
          TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
          (
            (1 - alphaX) * (1 - alphaY) * (- dnorm(zetaY[bY]) * pnorm((zetaX[bX + 1] - rhoAbn1 * zetaY[bY]) / sqrt(1 - rhoAbn1^2)) + 
                                             dnorm(zetaY[bY]) * pnorm((zetaX[bX] - rhoAbn1 * zetaY[bY]) / sqrt(1 - rhoAbn1^2))) + 
              (alphaX) * (1 - alphaY) * (- dnorm(zetaY[bY]) * pnorm((zetaX[bX + 1] - muX - rhoAbn3 * zetaY[bY]) / sqrt(1 - rhoAbn3^2)) + 
                                           dnorm(zetaY[bY]) * pnorm((zetaX[bX] - muX - rhoAbn3 * zetaY[bY]) / sqrt(1 - rhoAbn3^2))) + 
              (1 - alphaX) * (alphaY) * (- dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX + 1] - rhoAbn4 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn4^2)) + 
                                           dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX] - rhoAbn4 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn4^2))) + 
              (alphaX) * (alphaY) * (- dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX + 1] - muX - rhoAbn2 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn2^2)) + 
                                       dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX] - muX - rhoAbn2 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn2^2)))
          )
        # if (bY == 2){
        #   dev1 <- FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
        #     (
        #       - dnorm(zetaY[bY]) * pnorm((zetaX[bX + 1] - rhoNor * zetaY[bY]) / sqrt(1 - rhoNor^2)) + 
        #         dnorm(zetaY[bY]) * pnorm((zetaX[bX] - rhoNor * zetaY[bY]) / sqrt(1 - rhoNor^2))
        #     ) + 
        #     TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
        #                           (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
        #                           (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
        #                           alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        #     (
        #       (1 - alphaX) * (1 - alphaY) * (- dnorm(zetaY[bY]) * pnorm((zetaX[bX + 1] - rhoAbn1 * zetaY[bY]) / sqrt(1 - rhoAbn1^2)) + 
        #                                        dnorm(zetaY[bY]) * pnorm((zetaX[bX] - rhoAbn1 * zetaY[bY]) / sqrt(1 - rhoAbn1^2))) + 
        #         (alphaX) * (1 - alphaY) * (- dnorm(zetaY[bY]) * pnorm((zetaX[bX + 1] - muX - rhoAbn3 * zetaY[bY]) / sqrt(1 - rhoAbn3^2)) + 
        #                                      dnorm(zetaY[bY]) * pnorm((zetaX[bX] - muX - rhoAbn3 * zetaY[bY]) / sqrt(1 - rhoAbn3^2))) + 
        #         (1 - alphaX) * (alphaY) * (- dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX + 1] - rhoAbn4 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn4^2)) + 
        #                                      dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX] - rhoAbn4 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn4^2))) + 
        #         (alphaX) * (alphaY) * (- dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX + 1] - muX - rhoAbn2 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn2^2)) + 
        #                                  dnorm(zetaY[bY] - muY) * pnorm((zetaX[bX] - muX - rhoAbn2 * (zetaY[bY] - muY)) / sqrt(1 - rhoAbn2^2)))
        #     )
        #   dev2 <- (FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY] + delta), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) -
        #              FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))) / delta +
        #     (TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY] + delta), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) +
        #                               (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY] + delta), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) +
        #                               (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY] + delta), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) +
        #                               alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY] + delta), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))
        #      - TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) +
        #                                 (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) +
        #                                 (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) +
        #                                 alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))) / delta
        #   print(c(dev1, dev2))
        # }
      }
      
      if (bY < nBins){
        grad[10 + bY] <- grad[10 + bY] + 
          FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
          (
            dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX + 1] - rhoNor * zetaY[bY + 1]) / sqrt(1 - rhoNor^2)) - 
              dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX] - rhoNor * zetaY[bY + 1]) / sqrt(1 - rhoNor^2))
          ) + 
          TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
                                (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
                                (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
                                alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
          (
            (1 - alphaX) * (1 - alphaY) * (dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX + 1] - rhoAbn1 * zetaY[bY + 1]) / sqrt(1 - rhoAbn1^2)) - 
                                             dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX] - rhoAbn1 * zetaY[bY + 1]) / sqrt(1 - rhoAbn1^2))) + 
              (alphaX) * (1 - alphaY) * (dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX + 1] - muX - rhoAbn3 * zetaY[bY + 1]) / sqrt(1 - rhoAbn3^2)) - 
                                           dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX] - muX - rhoAbn3 * zetaY[bY + 1]) / sqrt(1 - rhoAbn3^2))) + 
              (1 - alphaX) * (alphaY) * (dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX + 1] - rhoAbn4 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn4^2)) - 
                                           dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX] - rhoAbn4 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn4^2))) + 
              (alphaX) * (alphaY) * (dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX + 1] - muX - rhoAbn2 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn2^2)) - 
                                       dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX] - muX - rhoAbn2 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn2^2)))
          )
        # if (bY == 1){
        #   dev1 <- FPCounts[bX, bY] / (pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)) * 
        #     (
        #       dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX + 1] - rhoNor * zetaY[bY + 1]) / sqrt(1 - rhoNor^2)) - 
        #         dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX] - rhoNor * zetaY[bY + 1]) / sqrt(1 - rhoNor^2))
        #     ) + 
        #     TPCounts[bX, bY] / ((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
        #                           (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
        #                           (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
        #                           alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)) * 
        #     (
        #       (1 - alphaX) * (1 - alphaY) * (dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX + 1] - rhoAbn1 * zetaY[bY + 1]) / sqrt(1 - rhoAbn1^2)) - 
        #                                        dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX] - rhoAbn1 * zetaY[bY + 1]) / sqrt(1 - rhoAbn1^2))) + 
        #         (alphaX) * (1 - alphaY) * (dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX + 1] - muX - rhoAbn3 * zetaY[bY + 1]) / sqrt(1 - rhoAbn3^2)) - 
        #                                      dnorm(zetaY[bY + 1]) * pnorm((zetaX[bX] - muX - rhoAbn3 * zetaY[bY + 1]) / sqrt(1 - rhoAbn3^2))) + 
        #         (1 - alphaX) * (alphaY) * (dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX + 1] - rhoAbn4 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn4^2)) - 
        #                                      dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX] - rhoAbn4 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn4^2))) + 
        #         (alphaX) * (alphaY) * (dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX + 1] - muX - rhoAbn2 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn2^2)) - 
        #                                  dnorm(zetaY[bY + 1] - muY) * pnorm((zetaX[bX] - muX - rhoAbn2 * (zetaY[bY + 1] - muY)) / sqrt(1 - rhoAbn2^2)))
        #     )
        # 
        #   dev2 <- (FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1] + delta), sigma = sigmaNor)) -
        #              FPCounts[bX, bY] * log(pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor))) / delta +
        #     (TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1] + delta), mean = c(0, 0), sigma = sigmaAbn1) +
        #                               (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1] + delta), mean = c(muX, 0), sigma = sigmaAbn3) +
        #                               (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1] + delta), mean = c(0, muY), sigma = sigmaAbn4) +
        #                               alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1] + delta), mean = c(muX, muY), sigma = sigmaAbn2))
        #      - TPCounts[bX, bY] * log((1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) +
        #                                 (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) +
        #                                 (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) +
        #                                 alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2))) / delta
        #   print(c(dev1, dev2))
        # }
      }
    }
  }
  
  return(-grad)
}

Dmvnorm <- function(XY, mean = c(0, 0), sigma){
  if (all(!is.finite(XY))){
    return(0)
  }
  else{
    x <- XY[1]
    y <- XY[2]
    return(dmvnorm(c(x, y), mean = mean, sigma = sigma))
  }
}


