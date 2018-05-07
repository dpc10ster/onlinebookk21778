# PlotROC.R
plotROC <- function (mu, sigma, FPF, TPF)
{
  zeta <- seq(- mu - 2,+ mu + 2,0.01)
  fpf <- array(dim = length(zeta))
  tpf <- array(dim = length(zeta))
  for (i in 1:length(zeta)) {
    fpf[i] <- pnorm(-zeta[i])
    tpf[i] <- pnorm((mu -zeta[i])/sigma) 
  }
  curveData <- data.frame(FPF = c(1, fpf, 0), TPF = c(1, tpf, 0))
  pointsData <- data.frame(FPF = FPF, TPF = TPF)
  rocPlot <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
    geom_line(data = curveData, size = 2) + geom_point(data = pointsData, shape = 1, size = 5)
  print(rocPlot)
}

