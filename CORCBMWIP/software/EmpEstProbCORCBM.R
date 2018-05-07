EmpEstProbCORCBM <- function(muX, muY, alphaX, alphaY, rhoNor, rhoAbn2, zetaX, zetaY , FPCounts, TPCounts){
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
  FPPrEst <- array(dim = c(nBins, nBins))
  TPPrEst <- FPPrEst
  for (bX in 1:nBins){
    for (bY in 1:nBins){
      FPPrEst[bX, bY] <- pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), sigma = sigmaNor)
      TPPrEst[bX, bY] <- (1 - alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, 0), sigma = sigmaAbn1) + 
        (alphaX) * (1 - alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, 0), sigma = sigmaAbn3) + 
        (1 - alphaX) * (alphaY) * pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(0, muY), sigma = sigmaAbn4) + 
        alphaX * alphaY* pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn2)
    }
  }
  
  FPCountsArray <- as.vector(FPCounts)
  FPOrder <- order(FPCountsArray, decreasing = TRUE)
  FPCountsArray <- FPCountsArray[FPOrder]
  FPPrArray <- as.vector(FPPrEst)
  FPPrArray <- FPPrArray[FPOrder]
  FPGroup <- NULL
  FPEstGroup <- NULL
  tempFP <- 0
  tempFPPr <- 0
  
  TPCountsArray <- as.vector(TPCounts)
  TPOrder <- order(TPCountsArray, decreasing = TRUE)
  TPCountsArray <- TPCountsArray[TPOrder]
  TPPrArray <- as.vector(TPPrEst)
  TPPrArray <- TPPrArray[TPOrder]
  TPGroup <- NULL
  TPEstGroup <- NULL
  tempTP <- 0
  tempTPPr <- 0
  for (ib in 1:nBins^2){
    if (FPCountsArray[ib] >= 5){
      FPGroup <- c(FPGroup, FPCountsArray[ib])
      FPEstGroup <- c(FPEstGroup, FPPrArray[ib])
    }else{
      tempFP <- tempFP + FPCountsArray[ib]
      tempFPPr <- tempFPPr + FPPrArray[ib]
      if (tempFP >= 5){
        FPGroup <- c(FPGroup, tempFP)
        FPEstGroup <- c(FPEstGroup, tempFPPr)
        tempFP <- 0
        tempFPPr <- 0
      }else if (ib == nBins^2){
        FPGroup[length(FPGroup)] <- FPGroup[length(FPGroup)] + tempFP
        FPEstGroup[length(FPEstGroup)] <- FPEstGroup[length(FPEstGroup)] + (1 - sum(FPEstGroup))
      }
    }
    
    if (TPCountsArray[ib] >= 5){
      TPGroup <- c(TPGroup, TPCountsArray[ib])
      TPEstGroup <- c(TPEstGroup, TPPrArray[ib])
    }else{
      tempTP <- tempTP + TPCountsArray[ib]
      tempTPPr <- tempTPPr + TPPrArray[ib]
      if (tempTP >= 5){
        TPGroup <- c(TPGroup, tempTP)
        TPEstGroup <- c(TPEstGroup, tempTPPr)
        tempTP <- 0
        tempTPPr <- 0
      }else if (ib == nBins^2){
        TPGroup[length(TPGroup)] <- TPGroup[length(TPGroup)] + tempTP
        TPEstGroup[length(TPEstGroup)] <- TPEstGroup[length(TPEstGroup)] + (1 - sum(TPEstGroup))
      }
    }
  }
  
  
  FPEmpGroup <- FPGroup / sum(FPGroup)
  TPEmpGroup <- TPGroup / sum(TPGroup)
  
  chiStat <- sum(FPGroup) * sum((FPEmpGroup - FPEstGroup)^2 / FPEstGroup) + sum(TPGroup) * sum((TPEmpGroup - TPEstGroup)^2 / TPEstGroup)
  df <- (length(FPEmpGroup) - 1) + (length(TPEmpGroup) - 1) - (6 + 2 * (nBins - 1))
  pVal <- 1 - pchisq(chiStat, df)
  
  # plot(FPPrEmp, FPPrEst)
  # lines(c(0, max(c(FPPrEmp, FPPrEst))), c(0, max(c(FPPrEmp, FPPrEst))))
  # plot(TPPrEmp, TPPrEst)
  # lines(c(0, max(c(TPPrEmp, TPPrEst))), c(0, max(c(TPPrEmp, TPPrEst))))
  return(list(FPEmpGroup = FPEmpGroup,
              FPEstGroup = FPEstGroup, 
              TPEmpGroup = TPEmpGroup,
              TPEstGroup = TPEstGroup,
              chiStat = chiStat,
              df = df,
              pVal = pVal))
  
}