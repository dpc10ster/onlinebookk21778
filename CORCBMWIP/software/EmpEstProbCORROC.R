EmpEstProbCORROC <- function(aX, bX, aY, bY, rho1, rho2, zetaX, zetaY , FPCounts, TPCounts){
  zetaX <- c(-Inf, zetaX, Inf)
  zetaY <- c(-Inf, zetaY, Inf)
  
  muX <- aX / bX
  muY <- aY / bY
  stdX <- 1 / bX
  stdY <- 1 / bY
  
  sigmaNor <- rbind(c(1, rho1), c(rho1, 1))
  sigmaAbn <- rbind(c(stdX^2, rho2 * stdX * stdY), c(rho2 * stdX * stdY, stdY^2))
  
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
      TPPrEst[bInx] <- pmvnorm(c(zetaX[bX], zetaY[bY]), c(zetaX[bX + 1], zetaY[bY + 1]), mean = c(muX, muY), sigma = sigmaAbn)
      TPNum[bInx] <- TPCounts[bX, bY]
      bInx <- bInx + 1
    }
  }
  FPPrEmp <- FPNum / sum(FPNum)
  TPPrEmp <- TPNum / sum(TPNum)
  
  chiStat <- sum(FPCounts) * sum((FPPrEmp - FPPrEst)^2 / FPPrEst) + sum(TPCounts) * sum((TPPrEmp - TPPrEst)^2 / TPPrEst)
  df <- (nBins^2 - 1) * 2 - (6 + 2 * (nBins - 1))
  pVal <- 1 - pchisq(chiStat, df)
  
  # plot(FPPrEmp, FPPrEst)
  # lines(c(0, max(c(FPPrEmp, FPPrEst))), c(0, max(c(FPPrEmp, FPPrEst))))
  # plot(TPPrEmp, TPPrEst)
  # lines(c(0, max(c(TPPrEmp, TPPrEst))), c(0, max(c(TPPrEmp, TPPrEst))))
  return(list(FPPrEmp = FPPrEmp,
              FPPrEst = FPPrEst, 
              TPPrEmp = TPPrEmp,
              TPPrEst = TPPrEst,
              chiStat = chiStat,
              df = df,
              pVal = pVal))
  
}