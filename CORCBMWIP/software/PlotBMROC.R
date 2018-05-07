PlotBMROC <- function(retCORROC, fpf, tpf){
  plotZeta <- seq(-20, 20, by = 0.01)
  
  if (retCORROC$done == 1){
    parms <- retCORROC$parms
    aX <- parms[1]
    bX <- parms[2]
    aY <- parms[3]
    bY <- parms[4]
    FPFX <- 1 - pnorm(plotZeta)
    FPFY <- 1 - pnorm(plotZeta)
    TPFX <- pnorm(aX - bX * plotZeta)
    TPFY <- pnorm(aY - bY * plotZeta)
    
    plotCORROC1 <- data.frame(FPF = FPFX, TPF = TPFX)
    plotCORROC2 <- data.frame(FPF = FPFY, TPF = TPFY)
    
    plotCurve1 <- plotCORROC1
    plotCurve2 <- plotCORROC2
  }
  fpf <- fpf[ , -dim(fpf)[2]]
  tpf <- tpf[ , -dim(tpf)[2]]
  
  
  plotOp1 <- data.frame(FPF = fpf[1, ], TPF = tpf[1, ], Modality = "1")
  FPF <- plotOp1$FPF
  TPF <- plotOp1$TPF
  
  fitPlot1 <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF), data = plotCurve1, size = 1) + 
    geom_point(mapping = aes(x = FPF, y = TPF), data = plotOp1, size = 3) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  
  return(list(fitPlot1))
}