PlotRSMBM <- function(a, b, mu, lambda, nu, nLesDistr, RowString){
  plotZeta <- seq(-20, 20, by = 0.01)
  lambdaP <- lambda / mu
  if (abs(nu * mu) <= 1e-6 ) nuP <- 1e-6 else nuP <- (1-exp(-nu * mu))
  
  FPFBM <- 1 - pnorm(plotZeta)
  TPFBM <- 1 - pnorm(plotZeta, mean = a/b, sd = 1/b)
  plotBM <- data.frame(FPF = FPFBM, TPF = TPFBM, Model = "BM")
  
  FPFSM <- sapply(plotZeta, xROC, lambdaP = lambdaP)
  TPFSM <- sapply(plotZeta, yROC, mu = mu, lambdaP = lambdaP, 
                  nuP = nuP, lesionDistribution = nLesDistr)
  plotSM <- data.frame(FPF = FPFSM, TPF = TPFSM, Model = "RSM")
  
  plotCurve <- rbind(plotSM, plotBM)
  dashedSM <- data.frame(FPF = c(FPFSM[1], 1), TPF = c(TPFSM[1], 1), Model = "RSM")
  
  rowTitle <- paste0("RowString = ", RowString)
  rowTitle <- RowString
  compPlot <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF, color = Model), data = plotCurve, size = 1) + 
    geom_line(data = dashedSM, aes(x = FPF, y = TPF, color = Model), linetype = 3, size = 1) + scale_color_manual(values = c("red", "blue")) +
    theme(legend.direction = "horizontal", legend.position = c(0.7, 0.1)) + 
    labs(title = rowTitle) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  
  fpfDashedSM <- seq(FPFSM[1], 1, by = 0.01)
  tpfDashedSM <- TPFSM[1] + (fpfDashedSM - FPFSM[1]) * (1 - TPFSM[1]) / (1 - FPFSM[1])
  fpfSM <- c(rev(fpfDashedSM), FPFSM)
  tpfSM <- c(rev(tpfDashedSM), TPFSM)
  fpf <- seq(0, 1, by = 0.01)
  tpfDiff <- rep(NA, length(fpf))
  for (i in 1:length(fpf)){
    indxSM <- which.min(abs(fpf[i] - fpfSM))
    indxBM <- which.min(abs(fpf[i] - FPFBM))
    tpfDiff[i] <- tpfSM[indxSM] - TPFBM[indxBM]
  }
  if (max(abs(tpfDiff)) < 0.01){
    plotDiff <- data.frame(FPF = fpf, TPF = tpfDiff, Model = "Diff")
    diffPlot <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF, color = Model), data = plotDiff, size = 1) + 
      theme(legend.direction = "horizontal", legend.position = c(0.8, 0.8)) + 
      labs(title = rowTitle) +
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0))
    return(list(compPlot = compPlot,
                diffPlot = diffPlot))
  }else{
    return(list(compPlot = compPlot))
  }
}