CBMCurve <- function(mu, alpha){
  plotZeta <- seq(-20, 20, by = 0.01)
  FPF <- 1 - pnorm(plotZeta)
  FPF <- c(1, FPF, 0)
  TPF <- (1 - alpha) * (1 - pnorm(plotZeta)) + alpha * (1 - pnorm(plotZeta, mean = mu))
  TPF <- c(1, TPF, 0)
  plotCurve <- data.frame(FPF = FPF, TPF = TPF)
  return(plotCurve)
}