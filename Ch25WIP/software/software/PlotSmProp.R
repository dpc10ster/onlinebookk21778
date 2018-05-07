PlotRsmProp <- function(mu, lambdaP, nuP, nLesDistr, c, da, fpf, tpf, i, j){
  if (c < 0){
    plotZeta <- seq(da/(4 * c) * sqrt(1 + c^2), 20, by = 0.005)
  }else if(c > 0){
    plotZeta <- seq(-20, da/(4 * c) * sqrt(1 + c^2), by = 0.005)
  }else{
    plotZeta <- seq(-20, 20, by = 0.005)
  }
  FPFProp <- pnorm(-(1 - c) * plotZeta - da/2 * sqrt(1 + c^2)) + pnorm(-(1 - c) * plotZeta + da/(2*c) * sqrt(1 + c^2)) - ifelse(c >= 0, 1, 0)
  TPFProp <- pnorm(-(1 + c) * plotZeta + da/2 * sqrt(1 + c^2)) + pnorm(-(1 + c) * plotZeta + da/(2*c) * sqrt(1 + c^2)) - ifelse(c >= 0, 1, 0)
  plotProp <- data.frame(FPF = FPFProp, TPF = TPFProp, Model = "PROP")
  
  plotZeta <- seq(-20, 20, by = 0.005)
  FPFSm <- sapply(plotZeta, xROC, lambdaP = lambdaP)
  TPFSm <- sapply(plotZeta, yROC, mu = mu, lambdaP = lambdaP, 
                  nuP = nuP, pmfLesionDistribution = nLesDistr)
  plotSm <- data.frame(FPF = FPFSm, TPF = TPFSm, Model = "SM")
  
  plotCurve <- rbind(plotSm, plotProp)
  dashedSm <- data.frame(FPF = c(FPFSm[1], 1), TPF = c(TPFSm[1], 1), Model = "SM")
  plotOp <- data.frame(FPF = fpf, TPF = tpf)
  
  ij <- paste0("i = ", i, ", j = ", j)
  fitPlot <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF, linetype = Model), data = plotCurve, size = 1) + 
    geom_line(data = dashedSm, aes(x = FPF, y = TPF), linetype = 3, size = 1) + 
    geom_point(mapping = aes(x = FPF, y = TPF), data = plotOp, size = 3) +
    theme(legend.position = "none") + 
    labs(title = ij) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  return(fitPlot)
}