#' @export
UtilGoodnessOfFit <- function(binCounts, binProb, nParam){
  retComb1 <- CombBins(binCounts, binProb)
  retComb1 <- CombBins(retComb1$obs[c(2, 1), , drop = FALSE], retComb1$prob[c(2, 1), , drop = FALSE])
  obs1 <- retComb1$obs[c(2, 1), , drop = FALSE]; exp1 <- retComb1$prob[c(2, 1), , drop = FALSE] * rowSums(obs1)
  fpGoodness1 <- rbind(obs1[1, ], exp1[1, ])
  tpGoodness1 <- rbind(obs1[2, ], exp1[2, ])
  
  retComb2 <- CombBins(binCounts[c(2, 1), , drop = FALSE], binProb[c(2, 1), , drop = FALSE])
  retComb2 <- CombBins(retComb2$obs[c(2, 1), , drop = FALSE], retComb2$prob[c(2, 1), , drop = FALSE])
  obs2 <- retComb2$obs; exp2 <- retComb2$prob * rowSums(obs2)
  fpGoodness2 <- rbind(obs2[1, ], exp2[1, ])
  tpGoodness2 <- rbind(obs2[2, ], exp2[2, ])
  
  if (ncol(fpGoodness1) >= ncol(fpGoodness2)){
    fpGoodness <- fpGoodness1
    tpGoodness <- tpGoodness1
    nBinsComb <- ncol(fpGoodness1)
  }else{
    fpGoodness <- fpGoodness2
    tpGoodness <- tpGoodness2
    nBinsComb <- ncol(fpGoodness2)
  }
  
  
  if (nBinsComb > (nParam + 1)){
    chisq <- sum((fpGoodness[1, ] - fpGoodness[2, ])^2/fpGoodness[2, ]) + sum((tpGoodness[1, ] - tpGoodness[2, ])^2/tpGoodness[2, ])
    df <- nBinsComb - (nParam + 1)
    pVal <- pchisq(chisq, df, lower.tail = FALSE)
  }else{
    chisq <- NA
    pVal <- NA
    df <- NA
  }
  return(list(chisq = chisq,
              df = df,
              pVal = pVal
              ))
}

CombBins <- function(binned, prob){
  if (ncol(binned) > 1){
    less5Indx <- which(binned[1, ] < 5)
    if (length(less5Indx) > 0){
      less5 <- binned[, less5Indx, drop = FALSE]
      less5Prob <- prob[, less5Indx, drop = FALSE]
      binned <- binned[, -less5Indx, drop = FALSE]
      prob <- prob[, -less5Indx, drop = FALSE]
      while (sum(less5[1, ]) >= 5){
        if (sum(less5[1, ]) == 5){
          binned <- cbind(binned, rowSums(less5))
          prob <- cbind(prob, rowSums(less5Prob))
          less5 <- numeric(0)
          less5Prob <- numeric(0)
          break
        }else{
          maxIndx <- which.max(less5[1, ])
          tempBin <- less5[, maxIndx, drop = FALSE]
          tempProb <- less5Prob[, maxIndx, drop = FALSE]
          less5 <- less5[, -maxIndx, drop = FALSE]
          less5Prob <- less5Prob[, -maxIndx, drop = FALSE]
          while (tempBin[1, ] < 5){
            minIndx <- which.min(less5[1, ])
            tempBin <- tempBin + less5[, minIndx, drop = FALSE]
            tempProb <- tempProb + less5Prob[, minIndx, drop = FALSE]
            less5 <- less5[, -minIndx, drop = FALSE]
            less5Prob <- less5Prob[, -minIndx, drop = FALSE]
          }
          binned <- cbind(binned, tempBin)
          prob <- cbind(prob, tempProb)
        }
      }
      if (length(less5) > 0){
        minIndx <- which.min(binned[1, ])
        binned[ , minIndx] <- binned[ , minIndx] + rowSums(less5)
        prob[ , minIndx] <- prob[ , minIndx] + rowSums(less5Prob)
      }
    }
  }
  return(list(
    obs = binned, 
    prob = prob
  ))
}