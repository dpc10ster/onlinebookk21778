PlotCORCBMCORROC <- function(retCORCBM, retCORROC, fpf, tpf, ciPoint){
  plotZeta <- seq(-20, 20, by = 0.01)
  muX <- retCORCBM$muX
  muY <- retCORCBM$muY
  alphaX <- retCORCBM$alphaX
  alphaY <- retCORCBM$alphaY
  
  FPFX <- 1 - pnorm(plotZeta)
  TPFX <- (1 - alphaX) * (1 - pnorm(plotZeta)) + alphaX * (1 - pnorm(plotZeta, mean = muX))
  FPFY <- 1 - pnorm(plotZeta)
  TPFY <- (1 - alphaY) * (1 - pnorm(plotZeta)) + alphaY * (1 - pnorm(plotZeta, mean = muY))
  
  plotCORCBM1 <- data.frame(FPF = FPFX, TPF = TPFX, Model = "CORCBM", Modality = "1")
  plotCORCBM2 <- data.frame(FPF = FPFY, TPF = TPFY, Model = "CORCBM", Modality = "2")
  
  if (retCORROC$done == 1){
    parms <- retCORROC$parms
    aX <- parms[1]
    bX <- parms[2]
    aY <- parms[3]
    bY <- parms[4]
    TPFX <- pnorm(aX - bX * plotZeta)
    TPFY <- pnorm(aY - bY * plotZeta)
    
    plotCORROC1 <- data.frame(FPF = FPFX, TPF = TPFX, Model = "CORROC2", Modality = "1")
    plotCORROC2 <- data.frame(FPF = FPFY, TPF = TPFY, Model = "CORROC2", Modality = "2")
    
    plotCurve1 <- rbind(plotCORCBM1, plotCORROC1)
    plotCurve2 <- rbind(plotCORCBM2, plotCORROC2)
  }else{
    plotCurve1 <- plotCORCBM1
    plotCurve2 <- plotCORCBM2
  }
  fpf <- fpf[ , -dim(fpf)[2]]
  tpf <- tpf[ , -dim(tpf)[2]]
  
  
  plotOp1 <- data.frame(FPF = fpf[1, ], TPF = tpf[1, ], Modality = "1")
  FPF <- plotOp1$FPF
  TPF <- plotOp1$TPF
  ciX <- binom.confint(x = FPF * K1, n = K1, method = "exact")
  ciY <- binom.confint(x = TPF * K2, n = K2, method = "exact")
  ciXUpper <- ciX$upper
  ciXLower <- ciX$lower
  ciYUpper <- ciY$upper
  ciYLower <- ciY$lower
  
  fitPlot1 <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF, linetype = Model), data = plotCurve1, size = 1) + 
    geom_point(mapping = aes(x = FPF, y = TPF), data = plotOp1, size = 3) +
    theme(legend.direction = "horizontal", legend.position = c(0.7, 0.1)) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  for (i in c(ciPoint[1])){
    ciX <- data.frame(FPF = c(ciXUpper[i], ciXLower[i]), TPF = c(TPF[i], TPF[i]))
    ciY <- data.frame(FPF = c(FPF[i], FPF[i]), TPF = c(ciYUpper[i], ciYLower[i]))
    fitPlot1 <- fitPlot1 + geom_line(data = ciY, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = ciX, aes(x = FPF, y = TPF), color = "black")
    barRgt <- data.frame(FPF = c(ciXUpper[i], ciXUpper[i]), TPF = c(TPF[i] - 0.01, TPF[i] + 0.01))
    barLft <- data.frame(FPF = c(ciXLower[i], ciXLower[i]), TPF = c(TPF[i] - 0.01, TPF[i] + 0.01))
    barUp <- data.frame(FPF = c(FPF[i] - 0.01, FPF[i] + 0.01), TPF = c(ciYUpper[i], ciYUpper[i]))
    barBtm <- data.frame(FPF = c(FPF[i] - 0.01, FPF[i] + 0.01), TPF = c(ciYLower[i], ciYLower[i]))
    fitPlot1 <- fitPlot1 + geom_line(data = barRgt, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barLft, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barUp, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barBtm, aes(x = FPF, y = TPF), color = "black")
  }
  
  
  plotOp2 <- data.frame(FPF = fpf[2, ], TPF = tpf[2, ], Modality = "2")
  FPF <- plotOp2$FPF
  TPF <- plotOp2$TPF
  ciX <- binom.confint(x = FPF * K1, n = K1, method = "exact")
  ciY <- binom.confint(x = TPF * K2, n = K2, method = "exact")
  ciXUpper <- ciX$upper
  ciXLower <- ciX$lower
  ciYUpper <- ciY$upper
  ciYLower <- ciY$lower
  
  fitPlot2 <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF, linetype = Model), data = plotCurve2, size = 1) + 
    geom_point(mapping = aes(x = FPF, y = TPF), data = plotOp2, size = 3) +
    theme(legend.direction = "horizontal", legend.position = c(0.7, 0.1)) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  for (i in c(ciPoint[2])){
    ciX <- data.frame(FPF = c(ciXUpper[i], ciXLower[i]), TPF = c(TPF[i], TPF[i]))
    ciY <- data.frame(FPF = c(FPF[i], FPF[i]), TPF = c(ciYUpper[i], ciYLower[i]))
    fitPlot2 <- fitPlot2 + geom_line(data = ciY, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = ciX, aes(x = FPF, y = TPF), color = "black")
    barRgt <- data.frame(FPF = c(ciXUpper[i], ciXUpper[i]), TPF = c(TPF[i] - 0.01, TPF[i] + 0.01))
    barLft <- data.frame(FPF = c(ciXLower[i], ciXLower[i]), TPF = c(TPF[i] - 0.01, TPF[i] + 0.01))
    barUp <- data.frame(FPF = c(FPF[i] - 0.01, FPF[i] + 0.01), TPF = c(ciYUpper[i], ciYUpper[i]))
    barBtm <- data.frame(FPF = c(FPF[i] - 0.01, FPF[i] + 0.01), TPF = c(ciYLower[i], ciYLower[i]))
    fitPlot2 <- fitPlot2 + geom_line(data = barRgt, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barLft, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barUp, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barBtm, aes(x = FPF, y = TPF), color = "black")
  }
  return(list(fitPlot1,
              fitPlot2))
}