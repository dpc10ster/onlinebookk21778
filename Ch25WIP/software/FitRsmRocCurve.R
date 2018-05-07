FitRsmRocCurve <- function(simuData, selectedMod, selectedRdr, AUCCbm, zetaCbm, nLesDistr){
  fpBinned <- simuData$NL
  tpBinned <- simuData$LL
  I <- dim(fpBinned)[1]
  J <- dim(tpBinned)[2]
  
  K <- dim(fpBinned)[3]
  K2 <- dim(tpBinned)[3]
  K1 <- K - K2
  modalityID <- simuData$modalityID
  readerID <- simuData$readerID
  
  I <- length(selectedMod)
  if (all(is.character(selectedMod)))
    selectedMod <- which( modalityID %in% selectedMod )
  fpBinned <- fpBinned[selectedMod, , 1:K1, ]
  tpBinned <- tpBinned[selectedMod, , , ]
  dim(fpBinned) <- c(I, J, K1, 1)
  dim(tpBinned) <- c(I, J, K2, 1)
  modalityID <- modalityID[selectedMod]
  
  J <- length(selectedRdr)
  if (all(is.character(selectedRdr)))
    selectedRdr <- which( readerID %in% selectedRdr )
  fpBinned <- fpBinned[ , selectedRdr, , ]
  tpBinned <- tpBinned[ , selectedRdr, , ]
  dim(fpBinned) <- c(I, J, K1, 1)
  dim(tpBinned) <- c(I, J, K2, 1)
  readerID <- readerID[selectedRdr]
  
  ROCPoints <- NULL
  ROCDashes <- NULL
  ROCOpPoints <- NULL
  lambdaOpt <- data.frame(lambda = rep(NA, prod(I, J)), Modality = NA, Reader = NA) 
  lambdaPOpt <- data.frame(lambdaP = rep(NA, prod(I, J)), Modality = NA, Reader = NA) 
  muOpt <- data.frame(mu = rep(NA, prod(I, J)), Modality = NA, Reader = NA) 
  nuPOpt <- data.frame(nuP = rep(NA, prod(I, J)), Modality = NA, Reader = NA)
  nuOpt <- data.frame(nu = rep(NA, prod(I, J)), Modality = NA, Reader = NA)
  AUC <- data.frame(AUC = rep(NA, prod(I, J)), Modality = NA, Reader = NA) 
  gdnss <- data.frame(gdnss = rep(NA, prod(I, J)), Modality = NA, Reader = NA) 
  plotIndex <- 1
  avg <- FALSE
  for (i in 1:I){
    for (j in 1:J){
      scores <- sort(unique(c(fpBinned[i, j, , ], tpBinned[i, j, , ])))
      nBins <- length(scores)
      fpBinnedTable <- rep(NA, nBins)
      tpBinnedTable <- fpBinnedTable
      for (b in 1:nBins){
        fpBinnedTable[b] <- sum(fpBinned == scores[b])
        tpBinnedTable[b] <- sum(tpBinned == scores[b])
      }
      
      FPF <- cumsum(rev(fpBinnedTable)) / K1
      TPF <- cumsum(rev(tpBinnedTable)) / K2
      lambdaPIni <- -log(1 - rev(FPF)[which(rev(FPF) < 1)[1]]) #??
      if (lambdaPIni == 0){
        errMsg <- paste0("Too few false positives on non-diseased cases in modality ", selectedMod[i], ", and reader ", selectedRdr[j], ".")
        stop(errMsg)
      }
      
      # lower bound of lambdaP. It limits the x-coord of the end point of the curve to be
      # larger than that of the second last (from upper rhs) operating point.
      lambdaPLower <- -log(1 - sort(unique(FPF), decreasing = TRUE)[3])
      
      minNuP <- RJafrocEnv$minNuP
      maxNuP <- RJafrocEnv$maxNuP
      if (FPF[length(FPF) - 1] >= TPF[length(FPF) - 1]){
        nuPIni <- RJafrocEnv$minNuP
      }else{
        TPFDiff <- function(lambdaP, nuP, nLesDistr, TPF){
          diff <- yROC(-Inf, 0, lambdaP, nuP, nLesDistr) - TPF
          return(diff)
        }
        nuPIni <- uniroot(TPFDiff, interval = c(minNuP, maxNuP), 
                          lambdaP = lambdaPIni, nLesDistr = nLesDistr, TPF = TPF[length(FPF) - 1], 
                          tol = .Machine$double.eps^0.5)$root
      }
      lambdaPIniFwd <- ForwardLambda(lambdaPIni, lambdaPLower)
      parameters <- lambdaPIniFwd
      namesVector <- "lambdaPFwd"
      names(parameters) <- namesVector
      zetaIni <- zetaCbm
      logLikIni <- iLoopVaryNuZetasRoc(lambdaPIni, fb = fpBinnedTable, tb = tpBinnedTable, K1 = K1, K2 = K2, nLesDistr = nLesDistr, 
                                    AUCCbm = AUCCbm[i, j], nuPIni = nuPIni, zetaIni = zetaIni, iniCal = TRUE)
      if (DEBUG) cat("logLikIni = ", logLikIni, ", lambdaPIni = ", lambdaPIni, "\n")
      # this defines the outer loop where lambdaFwd is varied
      ret <- mle2(oLoopVaryLambdaRoc, start = as.list(parameters),
                  data = list(fb = fpBinnedTable, tb = tpBinnedTable, K1 = K1, K2 = K2, nLesDistr = nLesDistr, 
                              AUCCbm = AUCCbm[i, j], nuPIni = nuPIni, lambdaPLower = lambdaPLower, zetaIni = zetaIni, 
                              iniCal = FALSE, weighted = FALSE), 
                  method = "BFGS") # vary lambdaP
      lambdaP <- as.numeric(InverseLamdba(ret@coef[1], lambdaPLower))
      logLikFin <- ret@min
      if (DEBUG) cat("logLikFin = ", logLikFin, ", lambdaP = ", lambdaP, "\n")
      # Look for other parameters by running inner loop using the optimized lambdaP
      ret <- iLoopVaryNuZetasRoc(lambdaP, fb = fpBinnedTable, tb = tpBinnedTable, K1 = K1, K2 = K2, nLesDistr = nLesDistr, 
                              AUCCbm = AUCCbm[i, j], nuPIni = nuPIni, zetaIni = zetaIni, iniCal = FALSE)
      nuPLower <- ret$nuPLower
      nuPUpper <- ret$nuPUpper
      ret <- ret$mleRet
      nuP <- as.numeric(InverseNu(ret@coef[1], nuPLower = nuPLower, nuPUpper = nuPUpper))
      
      mu <- uniroot(DiffAucSmMinusAucCbm, c(minMu, maxMu), lambdaP = lambdaP, nuP = nuP,
                    nLesDistr = nLesDistr, AUCCbm = AUCCbm[i, j], tol = .Machine$double.eps^0.5)$root
      finalZetas <- InverseZetas(ret@coef[2:length(ret@coef)])
      
      lambdaPOpt[plotIndex, ] <- c(lambdaP, modalityID[i], readerID[j])
      lambdaOpt[plotIndex, ] <- c(lambdaP * mu, modalityID[i], readerID[j])
      
      nuPOpt[plotIndex, ] <- c(nuP, modalityID[i], readerID[j])
      nu <- as.numeric(-log(1 - nuP) / abs(mu))
      nuOpt[plotIndex, ] <- c(nu, modalityID[i], readerID[j])
      
      muOpt[plotIndex, ] <- c(mu, modalityID[i], readerID[j])
      if (all(c(fpBinnedTable, tpBinnedTable) > 5)){
        gdnssTmp <- GdnsFtROC(mu, lambdaP, nuP, nLesDistr, fpBinnedTable, tpBinnedTable, K1, K2, finalZetas)
      }else{
        #message("Some bins have too few cases (less than 5). Goodness of fit is not returned.")
        gdnssTmp <- NA
      }
      gdnss[plotIndex, ] <- c(gdnssTmp, modalityID[i], readerID[j])
      if (!avg){
        plotStep <- 0.005
        plotZeta <- seq(from = -20, to = 20, by = plotStep)
        plotFPF <- sapply(plotZeta, xROC, lambdaP = as.numeric(lambdaPOpt[plotIndex, 1]))
        plotTPF <- sapply(plotZeta, yROC, mu = as.numeric(muOpt[plotIndex, 1]), lambdaP = as.numeric(lambdaPOpt[plotIndex, 1]), 
                          nuP = as.numeric(nuPOpt[plotIndex, 1]), pmfLesionDistribution = nLesDistr)
        deltaFPF <- plotFPF[1:(length(plotFPF) - 1)] - plotFPF[2:length(plotFPF)]   
        AUCTmp <- integrate(ROCInt, 0, max(plotFPF), mu = as.numeric(muOpt[plotIndex, 1]), lambdaP = as.numeric(lambdaPOpt[plotIndex, 1]), 
                            nuP = as.numeric(nuPOpt[plotIndex, 1]), pmfLesionDistribution = nLesDistr)$value
        AUC[plotIndex, ] <- c(AUCTmp + (1 + yROC(-20, as.numeric(muOpt[plotIndex, 1]), as.numeric(lambdaPOpt[plotIndex, 1]), as.numeric(nuPOpt[plotIndex, 1]), nLesDistr)) 
                              * (1 - max(plotFPF)) / 2,
                              modalityID[i], readerID[j])
        ROCPoints <- rbind(ROCPoints, data.frame(FPF = plotFPF, TPF = plotTPF, Modality = i, Reader = j))
        ROCDashes <- rbind(ROCDashes, data.frame(FPF = c(plotFPF[1], 1), TPF = c(plotTPF[1], 1), Modality = i, Reader = j))
        ROCOpPoints <- rbind(ROCOpPoints, data.frame(FPF = FPF[1:(length(FPF) - 1)], TPF = TPF[1:(length(TPF) - 1)], Modality = i, Reader = j))  
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
    lambdaPOpt <- data.frame(lambdaP = mean(as.numeric(lambdaPOpt$lambdaP)), Modality = paste(modalityID, collapse = "-"), Reader = paste(readerID, collapse = "-"))
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
  
  plots <- list()
  for (j in 1:length(selectedRdr)){
    plotTmp <- ggplot(data = ROCPoints[ROCPoints$Reader == readerID[j],], aes(x = FPF, y = TPF, color = Reader)) + 
      geom_line(size = 1, aes(linetype = Modality)) + scale_linetype_manual(values = c(1:(I+1))[-3])
    plotTmp <- plotTmp + geom_line(data = ROCDashes[ROCDashes$Reader == readerID[j],], aes(x = FPF, y = TPF, color = Reader, group = Modality), linetype = 3, size = 1) + 
      geom_point(data = ROCOpPoints[ROCOpPoints$Reader == readerID[j],], aes(shape = Modality), size = 5)
    ciX <- binom.confint(x = FPF * K1, n = K1, method = "exact")
    ciY <- binom.confint(x = TPF * K2, n = K2, method = "exact")
    ciXUpper <- ciX$upper
    ciXLower <- ciX$lower
    ciYUpper <- ciY$upper
    ciYLower <- ciY$lower
    for (i in c(1, (length(ROCOpPoints$FPF)))){
      ciX <- data.frame(FPF = c(ciXUpper[i], ciXLower[i]), TPF = c(TPF[i], TPF[i]))
      ciY <- data.frame(FPF = c(FPF[i], FPF[i]), TPF = c(ciYUpper[i], ciYLower[i]))
      plotTmp <- plotTmp + geom_line(data = ciY, aes(x = FPF, y = TPF), color = "black") + 
        geom_line(data = ciX, aes(x = FPF, y = TPF), color = "black")
      barRgt <- data.frame(FPF = c(ciXUpper[i], ciXUpper[i]), TPF = c(TPF[i] - 0.01, TPF[i] + 0.01))
      barLft <- data.frame(FPF = c(ciXLower[i], ciXLower[i]), TPF = c(TPF[i] - 0.01, TPF[i] + 0.01))
      barUp <- data.frame(FPF = c(FPF[i] - 0.01, FPF[i] + 0.01), TPF = c(ciYUpper[i], ciYUpper[i]))
      barBtm <- data.frame(FPF = c(FPF[i] - 0.01, FPF[i] + 0.01), TPF = c(ciYLower[i], ciYLower[i]))
      plotTmp <- plotTmp + geom_line(data = barRgt, aes(x = FPF, y = TPF), color = "black") + 
        geom_line(data = barLft, aes(x = FPF, y = TPF), color = "black") + 
        geom_line(data = barUp, aes(x = FPF, y = TPF), color = "black") + 
        geom_line(data = barBtm, aes(x = FPF, y = TPF), color = "black")
    }
    plots[[j]] <- plotTmp + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
      theme(legend.position = "none", legend.direction = "horizontal", legend.justification = c(1, 0))
  }
  
  return(list(
    plots = plots,
    mu = muOpt,
    lambda = lambdaOpt,
    lambdaP = lambdaPOpt,
    nu = nuOpt,
    nuP = nuPOpt,
    AUC = AUC, 
    zetas = finalZetas, 
    gdnss = gdnss
  ))
}

GdnsFtROC <- function(mu, lambdaP, nuP, nLesDistr, fpBinnedTable, tpBinnedTable, K1, K2, finalZetas){
  zeta <- finalZetas
  FPF <- xROC(zeta, lambdaP)
  TPF <- yROC(zeta, mu, lambdaP, nuP, nLesDistr)
  FPF <- c(1, FPF) - c(FPF, 0)
  TPF <- c(1, TPF) - c(TPF, 0)
  numFP <- FPF * K1
  numTP <- TPF * K2
  chiStat <- sum((numFP - (fpBinnedTable))^2/numFP, (numTP - (tpBinnedTable))^2/numTP)
  dgrfrdm <- length(fpBinnedTable) - 2
  pval <- 1 - pchisq(chiStat, dgrfrdm)
  return(pval)
}