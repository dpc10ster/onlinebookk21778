FitURSM <- function(simuData, selectedMod, selectedRdr, nLesDistr){ 
  maxNuP <- RJafrocEnv$maxNuP
  minNuP <- RJafrocEnv$minNuP
  modalityID <- simuData$modalityID
  readerID <- simuData$readerID
  NL <- simuData$NL
  LL <- simuData$LL
  maxNumMarks <- dim(simuData$NL)[4]
  lesionNum <- simuData$lesionNum
  I <- length(selectedMod)
  J <- dim(NL)[2]
  K <- dim(NL)[3]
  K2 <- dim(LL)[3]
  K1 <- K - K2  
  
  if (all(is.character(selectedMod)))
    selectedMod <- which( modalityID %in% selectedMod )
  NL <- NL[selectedMod, , , ]
  LL <- LL[selectedMod, , , ]
  dim(NL) <- c(I, J, K, maxNumMarks)
  dim(LL) <- c(I, J, K2, max(lesionNum))
  modalityID <- modalityID[selectedMod]
  J <- length(selectedRdr)
  if (all(is.character(selectedRdr)))
    selectedRdr <- which( readerID %in% selectedRdr )
  NL <- NL[ , selectedRdr, , ]
  LL <- LL[ , selectedRdr, , ]
  dim(NL) <- c(I, J, K, maxNumMarks)
  dim(LL) <- c(I, J, K2, max(lesionNum))
  readerID <- readerID[selectedRdr]
  
  fpBinned <- apply(NL, c(1, 2, 3), max)
  fpBinned <- fpBinned[ , , 1:K1]
  tpBinned <- apply(LL, c(1, 2, 3), max)
  dim(fpBinned) <- c(I, J, K1)
  dim(tpBinned) <- c(I, J, K2)
   
  ROCPoints <- NULL
  ROCDashes <- NULL
  ROCOpPoints <- NULL
  lambdaPOpt <- data.frame(lambdaP = rep(NA, prod(I, J)), Modality = NA, Reader = NA) 
  muOpt <- data.frame(mu = rep(NA, prod(I, J)), Modality = NA, Reader = NA) 
  nuPOpt <- data.frame(nuP = rep(NA, prod(I, J)), Modality = NA, Reader = NA)
  AUC <- data.frame(AUC = rep(NA, prod(I, J)), Modality = NA, Reader = NA)  
  avg <- FALSE
  plotIndex <- 1
  for (i in 1:I){
    for (j in 1:J){
      scores <- sort(unique(c(fpBinned[i, j, ], tpBinned[i, j, ])))
      nBins <- length(scores)
      fpBinnedTable <- rep(NA, nBins)
      tpBinnedTable <- fpBinnedTable
      for (b in 1:nBins){
        fpBinnedTable[b] <- sum(fpBinned == scores[b])
        tpBinnedTable[b] <- sum(tpBinned == scores[b])
      }
      FPF <- cumsum(rev(fpBinnedTable)) / K1
      TPF <- cumsum(rev(tpBinnedTable)) / K2
      FPF <- FPF[-length(FPF)]
      TPF <- TPF[-length(TPF)]
      iniPntIndx <- which(TPF > FPF)
      iniPntIndx <- iniPntIndx[length(iniPntIndx)]
      lambdaP <- -log(1 - FPF[iniPntIndx]) #??
      if (lambdaP == 0){
        errMsg <- paste0("Too few false positives on non-diseased cases in modality ", trts[i], ", and reader ", rdrs[j], ".")
        stop(errMsg)
      }
      
      nuP <- uniroot(TPFDiff, interval = c(minNuP, maxNuP), 
                        lambdaP = lambdaP, nLesDistr = nLesDistr, TPF = TPF[iniPntIndx], 
                        tol = .Machine$double.eps^0.5)$root
      
      nuP <- min(nuP, maxNuP - 2 * .Machine$double.eps^0.5)
      nuP <- max(nuP, minNuP + 2 * .Machine$double.eps^0.5)
      sumNLFit <- rev(cumsum(rev(fpBinnedTable)))      
      
      zetaIni <- rep(NA, length(fpBinnedTable) - 1)
      revFPF <- rev(cumsum(rev(fpBinnedTable)) / K1)[-1] # gets rid of the 1, leaving the highest non-trivial FPF value

      # following follows from expression for FPF and relation between erf and Phi functions
      phiZeta <- (log(1 - revFPF)/lambdaP + 1) # xz was right
      for (z in 2:length(phiZeta)){
        phiZeta[z] <- max(phiZeta[z - 1] + 0.01, phiZeta[z])
      }
      if (any(phiZeta < 0.01)) {
        # dpc!! if this error occurs, the trial lambda is too small
        for (z in 1:length(phiZeta)){
          if (phiZeta[z] < 0.01){
            # impossible to found Gaussian quantile
            if (z == 1){
              # assign a small probability
              phiZeta[z] <- 0.01
            }else{
              # assign a small probability and larger than the previous cutoff
              phiZeta[z] <- phiZeta[z - 1] + 0.01
            }
          }else{
            # make sure every cutoff is larger than the previous one
            phiZeta[z] <- max(phiZeta[z - 1] + 0.01, phiZeta[z])
          }
        }
      }

      if (any(phiZeta > 0.99)){
        for (z in length(phiZeta):1){
          if (phiZeta[z] > 0.99){
            if (z == length(phiZeta)){
              phiZeta[z] <- 0.99
            }else{
              phiZeta[z] <- phiZeta[z + 1] - 0.01
            }
          }else{
            phiZeta[z] <- min(phiZeta[z + 1] - 0.01, phiZeta[z])
          }
        }
      }
      zetaIni <- qnorm(phiZeta)
      # cbmRet <- FitCbmRoc(simuData, selectedMod, selectedRdr)
      # zetaIni <- cbmRet$cutoffs
     zetaIniFwd <- ForwardZetas(zetaIni)
      
      lambdaPIniFwd <- ForwardLambda(lambdaP, minLambdaP)
      nuPIniFwd <- ForwardNu(nuP, minNuP, maxNuP)
      parameters <- c(lambdaPIniFwd, nuPIniFwd, zetaIniFwd) 
      namesVector <- NULL
      for (t in 1:length(zetaIniFwd))
        namesVector <- c(namesVector, paste0("zetaFwd", t))
      names(parameters) <- c("lambdaPFwd", "nuPFwd", namesVector)
      NllROCNew <- addArguments(NllROC, length(zetaIniFwd))
      ret <- mle2(NllROCNew, start = as.list(parameters),
                  data = list(sumNLFit = sumNLFit, fb = fpBinnedTable, tb = tpBinnedTable, K1 = K1, K2 = K2, nLesDistr = nLesDistr), 
                  method = "BFGS", control = list(maxit = 200))
      parametersNor <- as.vector(InverseTransform(ret@coef))
      lambdaPOpt[plotIndex, ] <- c(InverseLamdba(ret@coef[1], minLambdaP), modalityID[i], readerID[j])
      nuPOp <- InverseNu(ret@coef[2], minNuP, maxNuP)
      # zetaROC <- sumNLFit
      # for (t in 1:length(sumNLFit)) 
      #   zetaROC[t] <- optimize(ChisqZetaROC, c(minZeta, maxZeta), lambdaP = as.numeric(lambdaPOpt[plotIndex, 1]), numNL = sumNLFit[t], K1 = K1, tol = 1e-15)$minimum
      zetaROC <- InverseZetas(ret@coef[3:length(ret@coef)])
      muOpt[plotIndex, ] <- c(optimize(ChisqMuROC, c(minMu, maxMu), 
                                    zeta = zetaROC, fb = fpBinnedTable, tb = tpBinnedTable, K1 = K1, K2 = K2, lambdaP = as.numeric(lambdaPOpt[plotIndex, 1]), nuP = nuPOp, nLesDistr)$minimum,
                           modalityID[i], readerID[j])
      nuPOpt[plotIndex, ] <- c(nuPOp, modalityID[i], readerID[j])
      if (!avg){
        plotStep <- 0.005
        plotZeta <- seq(from = -20, to = 20, by = plotStep)
        plotFPF <- sapply(plotZeta, xROC, lambdaP = as.numeric(lambdaPOpt[plotIndex, 1]))
        plotTPF <- sapply(plotZeta, yROC, mu = as.numeric(muOpt[plotIndex, 1]), lambdaP = as.numeric(lambdaPOpt[plotIndex, 1]), nuP = as.numeric(nuPOpt[plotIndex, 1]), nLesDistr)
        deltaFPF <- plotFPF[1:(length(plotFPF) - 1)] - plotFPF[2:length(plotFPF)]        
        AUC[plotIndex, ] <- c(sum((plotTPF[1:(length(plotTPF) - 1)] + plotTPF[2:length(plotTPF)]) * deltaFPF / 2) + (plotTPF[1] + 1) * (1 - plotFPF[1]) / 2,
                           modalityID[i], readerID[j])
        ROCPoints <- rbind(ROCPoints, data.frame(FPF = plotFPF, TPF = plotTPF, Modality = i, Reader = j))
        ROCDashes <- rbind(ROCDashes, data.frame(FPF = c(plotFPF[1], 1), TPF = c(plotTPF[1], 1), Modality = i, Reader = j))
        ROCOpPoints <- rbind(ROCOpPoints, data.frame(FPF = FPF, TPF = TPF, Modality = i, Reader = j))  
      }
      plotIndex <- plotIndex + 1
    }
  }
  
  if (!avg){
    class <- paste("M-", modalityID[ROCPoints$Modality],"\n", "R-", readerID[ROCPoints$Reader], sep = "")
    ROCPoints <- data.frame(FPF = ROCPoints$FPF, TPF = ROCPoints$TPF, class = class, Modality = modalityID[ROCPoints$Modality], Reader = readerID[ROCPoints$Reader], type = "individual")
    class <- paste("M-", modalityID[ROCDashes$Modality],"\n", "R-", readerID[ROCDashes$Reader], sep = "")
    ROCDashes <- data.frame(FPF = ROCDashes$FPF, TPF = ROCDashes$TPF, class = class, Modality = modalityID[ROCDashes$Modality], Reader = readerID[ROCDashes$Reader], type = "individual")
    class <- paste("M-", modalityID[ROCOpPoints$Modality],"\n", "R-", readerID[ROCOpPoints$Reader], sep = "")
    ROCOpPoints <- data.frame(FPF = ROCOpPoints$FPF, TPF = ROCOpPoints$TPF, class = class, Modality = modalityID[ROCOpPoints$Modality], Reader = readerID[ROCOpPoints$Reader], type = "individual")
  }else{
    lambdaOpt <- data.frame(lambda = mean(as.numeric(lambdaOpt$lambda)), Modality = paste(modalityID, collapse = "-"), Reader = paste(readerID, collapse = "-"))
    lambdaPOpt <- data.frame(lambda = mean(as.numeric(lambdaPOpt$lambdaP)), Modality = paste(modalityID, collapse = "-"), Reader = paste(readerID, collapse = "-"))
    muOpt <- data.frame(mu = mean(as.numeric(muOpt$mu)), Modality = paste(modalityID, collapse = "-"), Reader = paste(readerID, collapse = "-")) 
    nuPOpt <- data.frame(nuP = mean(as.numeric(nuPOpt$nuP)), Modality = paste(modalityID, collapse = "-"), Reader = paste(readerID, collapse = "-")) 
    nuOpt <- data.frame(nu = mean(as.numeric(nuOpt$nu)), Modality = paste(modalityID, collapse = "-"), Reader = paste(readerID, collapse = "-")) 
    plotStep <- 0.005
    plotZeta <- seq(from = -20, to = 20, by = plotStep)
    plotFPF <- sapply(plotZeta, xROC, lambdaP = as.numeric(lambdaPOpt[1, 1]))
    plotTPF <- sapply(plotZeta, yROC, mu = as.numeric(muOpt[1, 1]), 
                      lambdaP = as.numeric(lambdaPOpt[1, 1]), nuP = as.numeric(nuPOpt[1, 1]), nLesDistr = nLesDistr)
    deltaFPF <- plotFPF[1:(length(plotFPF) - 1)] - plotFPF[2:length(plotFPF)]
    AUC <- data.frame(AUC = sum((plotTPF[1:(length(plotTPF) - 1)] + plotTPF[2:length(plotTPF)]) * deltaFPF / 2) + (plotTPF[1] + 1) * (1 - plotFPF[1]) / 2,
                      Modality = paste(modalityID, collapse = "-"), Reader = paste(readerID, collapse = "-")) 
    class <- paste(paste("M-"), paste(modalityID, collapse = " "),"\n", paste("R-"), paste(readerID, collapse = " "), sep = "")
    ROCPoints <- data.frame(FPF = plotFPF, TPF = plotTPF, class = class, type = "averaged")
    ROCDashes <- data.frame(FPF = c(plotFPF[1], 1), TPF = c(plotTPF[1], 1), class = class, type = "averaged")  
  }
  lambdaPOpt$lambdaP <- as.numeric(lambdaPOpt$lambdaP)
  muOpt$mu <- as.numeric(muOpt$mu)
  nuPOpt$nuP <- as.numeric(nuPOpt$nuP)
  AUC$AUC <- as.numeric(AUC$AUC)
  plots <- list()
  for (j in 1:length(selectedRdr)){
    plotTmp <- ggplot(data = ROCPoints[ROCPoints$Reader == readerID[j],], aes(x = FPF, y = TPF, color = Reader)) + 
      geom_line(size = 1, aes(linetype = Modality)) + scale_linetype_manual(values = c(1:(I+1))[-3])
    plotTmp <- plotTmp + geom_line(data = ROCDashes[ROCDashes$Reader == readerID[j],], aes(x = FPF, y = TPF, color = Reader, group = Modality), linetype = 3, size = 1) + 
      geom_point(data = ROCOpPoints[ROCOpPoints$Reader == readerID[j],], aes(shape = Modality), size = 5) + 
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0))
    plots[[j]] <- list()
    plots[[j]] <- plotTmp
  }
  return(list(
    plots = plots, 
    lambdaP = lambdaPOpt,
    mu = muOpt,
    nuP = nuPOpt,
    AUC = AUC)
    ) 
}


ChisqZetaROC <- function(zeta, lambdaP, numNL, K1){
  prob <- xROC(zeta, lambdaP)
  if (any(prob < 1e-15)) prob[which(prob < 1e-15)] <- 1e-15
  if (any(prob > (1 - 1e-15))) prob[which(prob > (1 - 1e-15))] <- 1 - 1e-15
  rval <- (numNL - prob * K1)^2 / (prob * K1)
  prob <- 1 - prob  
  rval <- rval + (K1 - numNL - prob * K1)^2 / (prob * K1)
  return(rval)
}

ChisqMuROC <- function(mu, zeta, fb, tb, K1, K2, lambdaP, nuP, nLesDistr){
  FPF <- xROC(zeta, lambdaP)
  TPF <- yROC(zeta, mu, lambdaP, nuP, nLesDistr)
  FPFTerms <- c(1, FPF) - c(FPF, 0)
  TPFTerms <- c(1, TPF) - c(TPF, 0)
  FPFTerms[ FPFTerms < 1e-15 ] <-  1e-15
  TPFTerms[ TPFTerms < 1e-15 ] <-  1e-15
  FPFTerms <- FPFTerms * K1
  TPFTerms <- TPFTerms * K2
  rval <- sum((FPFTerms - fb)^2 / FPFTerms, (TPFTerms - tb)^2 / TPFTerms)
  return(
    rval
  )
}

NllROC <- function(lambdaPFwd, nuPFwd, sumNLFit, fb, tb, K1, K2, nLesDistr){  
  lambdaP <- InverseLamdba(lambdaPFwd, minLambdaP)
  nuP <- InverseNu(nuPFwd, minNuP, maxNuP)
  allParameters <- names(formals())
  zetaPos <- grep("zetaFwd", allParameters)
  zeta <- unlist(mget(allParameters[zetaPos]))
  zeta <- InverseZetas(zeta)
  
  mu <- optimize(ChisqMuROC, c(minMu, maxMu), 
                 zeta = zeta, fb = fb, tb = tb, K1 = K1, K2 = K2, lambdaP = lambdaP, nuP = nuP, nLesDistr = nLesDistr)$minimum
  FPF <- xROC(zeta, lambdaP)
  TPF <- yROC(zeta, mu, lambdaP, nuP, nLesDistr)
  FPFTerms <- c(1, FPF) - c(FPF, 0)
  TPFTerms <- c(1, TPF) - c(TPF, 0)
  FPFTerms[ FPFTerms < 1e-15 ] <-  1e-15
  TPFTerms[ TPFTerms < 1e-15 ] <-  1e-15
  L <- sum(log(c(FPFTerms, TPFTerms)) * c(fb, tb)) 
  return (-L)
}

TPFDiff <- function(lambdaP, nuP, nLesDistr, TPF){
  diff <- yROC(-Inf, 0, lambdaP, nuP, nLesDistr) - TPF
  return(diff)
}
