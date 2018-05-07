PlotCBMROC <- function(mu, alpha){
  plotZeta <- seq(-20, 20, by = 0.01)
  FPF <- 1 - pnorm(plotZeta)
  TPF <- (1 - alpha) * (1 - pnorm(plotZeta)) + alpha * (1 - pnorm(plotZeta, mean = mu))
  
  plotCORCBM <- data.frame(FPF = FPF, TPF = TPF)
  cbmROCCurve <- ggplot() + geom_line(mapping = aes(x = FPF, y = TPF), data = plotCORCBM, size = 1) 
  vhvOp <- data.frame(x = c(0), y = c(0.75))
  cbmROCCurve <- cbmROCCurve + geom_point(data = vhvOp)
  return(cbmROCCurve)
}