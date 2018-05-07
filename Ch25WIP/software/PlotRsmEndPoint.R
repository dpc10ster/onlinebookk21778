PlotRsmEndPoint <- function(mu, lambdaP, nuP, nLesDistr){
  plotZeta <- seq(-30, -6, by = 0.005)
  FPFRsm <- sapply(plotZeta, xROC, lambdaP = lambdaP)
  TPFRsm <- sapply(plotZeta, yROC, mu = mu, lambdaP = lambdaP, 
                  nuP = nuP, pmfLesionDistribution = nLesDistr)
  plotRsm <- data.frame(FPF = FPFRsm, TPF = TPFRsm)
  
  k <- (1 - TPFRsm[1]) / (1 - FPFRsm[1])
  deltaFPF <- 1e-9
  dashedRsm <- data.frame(FPF = c(FPFRsm[1], FPFRsm[1] + deltaFPF), TPF = c(TPFRsm[1], TPFRsm[1] + k*deltaFPF))
  fitPlot <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF), data = plotRsm, size = 1) + 
    geom_line(data = dashedRsm, aes(x = FPF, y = TPF), linetype = 3, size = 1) + scale_color_manual(values = c("red")) +
    theme(legend.position = "none") + 
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  #cat("Num SM AUC =", trapz(c(rev(FPFRsm), 1), c(rev(TPFRsm), 1)), "Num PROP AUC =", trapz(rev(FPFProp), rev(TPFProp)), "\n")
  return(fitPlot)
}