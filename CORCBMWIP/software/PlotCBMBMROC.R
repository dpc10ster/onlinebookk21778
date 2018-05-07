PlotCBMBMROC <- function(retCORCBM, retCORROC, fpf, tpf, ciPoint){
  plotZeta <- seq(-20, 20, by = 0.01)
  muX <- retCORCBM$muX
  muY <- retCORCBM$muY
  alphaX <- retCORCBM$alphaX
  alphaY <- retCORCBM$alphaY
  
  FPFX <- 1 - pnorm(plotZeta)
  TPFX <- (1 - alphaX) * (1 - pnorm(plotZeta)) + alphaX * (1 - pnorm(plotZeta, mean = muX))
  FPFY <- 1 - pnorm(plotZeta)
  TPFY <- (1 - alphaY) * (1 - pnorm(plotZeta)) + alphaY * (1 - pnorm(plotZeta, mean = muY))
  
  plotCORCBM1 <- data.frame(FPF = FPFY, TPF = TPFY, Model = "CBM", Modality = "1")
  
  if (retCORROC$done == 1){
    parms <- retCORROC$parms
    aX <- parms[1]
    bX <- parms[2]
    aY <- parms[3]
    bY <- parms[4]
    TPFX <- pnorm(aX - bX * plotZeta)
    TPFY <- pnorm(aY - bY * plotZeta)
    
    plotCORROC1 <- data.frame(FPF = FPFY, TPF = TPFY, Model = "Binormal", Modality = "1")
    
    plotCurve1 <- rbind(plotCORCBM1, plotCORROC1)
  }else{
    plotCurve1 <- plotCORCBM1
  }
  fpf <- fpf[ , -dim(fpf)[2]]
  tpf <- tpf[ , -dim(tpf)[2]]
  
  
  plotOp1 <- data.frame(FPF = fpf[2, ], TPF = tpf[2, ], Modality = "1")
  FPF <- plotOp1$FPF
  TPF <- plotOp1$TPF
  
  fitPlot1 <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF, linetype = Model), data = plotCurve1, size = 1) + 
    geom_point(mapping = aes(x = FPF, y = TPF), data = plotOp1, size = 3) +
    theme(legend.direction = "horizontal", legend.position = c(0.7, 0.1)) + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  
  return(fitPlot1)
}