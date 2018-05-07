BMCurve <- function(a, b){
  plotZeta <- seq(-20, 20, by = 0.01)
  FPF <- 1 - pnorm(plotZeta)
  FPF <- c(1, FPF, 0)
  TPF <- pnorm(a - b * plotZeta)  
  TPF <- c(1, TPF, 0)
  plotCurve <- data.frame(FPF = FPF, TPF = TPF)
  return(plotCurve)
}