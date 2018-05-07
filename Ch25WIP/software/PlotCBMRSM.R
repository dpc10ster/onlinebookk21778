PlotCBMRSM <- function(muCbm, alpha, muSm, lambdaP, nuP, nLesDistr, fpf, tpf, i, j){
  plotZeta <- seq(-20, 20, by = 0.01)
  
  FPFCBM <- 1 - pnorm(plotZeta)
  TPFCBM <- (1 - alpha) * (1 - pnorm(plotZeta)) + alpha * (1 - pnorm(plotZeta, mean = muCbm))
  plotCBM <- data.frame(FPF = FPFCBM, TPF = TPFCBM, Model = "CBM")
  
  FPFSM <- sapply(plotZeta, xROC, lambdaP = lambdaP)
  TPFSM <- sapply(plotZeta, yROC, mu = muSm, lambdaP = lambdaP, 
                  nuP = nuP, pmfLesionDistribution = nLesDistr)
  plotSM <- data.frame(FPF = FPFSM, TPF = TPFSM, Model = "RSM")
  
  plotCurve <- rbind(plotSM, plotCBM)
  dashedSM <- data.frame(FPF = c(FPFSM[1], 1), TPF = c(TPFSM[1], 1), Model = "RSM")
  plotOp <- data.frame(FPF = fpf, TPF = tpf)
  
  ij <- paste0("i = ", i, ", j = ", j)
  fitPlot <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF, color = Model), data = plotCurve, size = 1) + 
    geom_line(data = dashedSM, aes(x = FPF, y = TPF, color = Model), linetype = 3, size = 1) + scale_color_manual(values = c("red", "blue")) +
    geom_point(mapping = aes(x = FPF, y = TPF), data = plotOp, size = 3) +
    theme(legend.direction = "horizontal", legend.position = c(0.7, 0.1)) + 
    labs(title = ij) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  return(fitPlot)
}