PlotRsmProp <- function(mu, lambdaP, nuP, nLesDistr, c1, da, fpf, tpf, i, j, K1, K2, errBar){
  if (c1 < 0){
    plotZeta <- seq(da/(4 * c1) * sqrt(1 + c1^2), mu + 4, by = 0.1)
  }else if(c1 > 0){
    plotZeta <- seq(-3, da/(4 * c1) * sqrt(1 + c1^2), by = 0.1)
  }else{
    plotZeta <- seq(-3, mu + 4, by = 0.1)
  }
  FPFProp <- pnorm(-(1 - c1) * plotZeta - da/2 * sqrt(1 + c1^2)) + pnorm(-(1 - c1) * plotZeta + da/(2*c1) * sqrt(1 + c1^2)) - ifelse(c1 >= 0, 1, 0)
  TPFProp <- pnorm(-(1 + c1) * plotZeta + da/2 * sqrt(1 + c1^2)) + pnorm(-(1 + c1) * plotZeta + da/(2*c1) * sqrt(1 + c1^2)) - ifelse(c1 >= 0, 1, 0)
  plotProp <- data.frame(FPF = c(FPFProp), TPF = c(TPFProp), Model = "PROP")
  prop11 <- data.frame(FPF = c(max(FPFProp), 1), TPF = c(max(TPFProp), 1), Model = "PROP")
  
  plotZeta <- seq(-3, mu + 4, by = 0.1)
  FPFSm <- sapply(plotZeta, xROC, lambdaP = lambdaP)
  TPFSm <- sapply(plotZeta, yROC, mu = mu, lambdaP = lambdaP, 
                  nuP = nuP, pmfLesionDistribution = nLesDistr)
  plotSm <- data.frame(FPF = FPFSm, TPF = TPFSm, Model = "URSM")
  
  plotCurve <- rbind(plotProp, plotSm)
  dashedSm <- data.frame(FPF = c(FPFSm[1], 1), TPF = c(TPFSm[1], 1), Model = "URSM")
  plotOp <- data.frame(FPF = fpf, TPF = tpf)
  
  ij <- paste0("i = ", i, ", j = ", j)
  fitPlot <- ggplot() + geom_line(data = prop11, aes(x = FPF, y = TPF, color = Model), size = 1) +
    geom_line(mapping = aes(x = FPF, y = TPF, color = Model), data = plotCurve, size = 1) + 
    geom_line(data = dashedSm, aes(x = FPF, y = TPF, color = Model), linetype = 3, size = 1) + 
    scale_color_manual(values = c("blue", "red")) +
    geom_point(mapping = aes(x = FPF, y = TPF), data = plotOp, size = 3) +
    theme(legend.position = "none") + 
    labs(title = ij) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  
  if (length(errBar) == 1 && errBar == "all") {
    ciIndx <- TRUE
  }else {
    ciIndx <- errBar
  }
  FPF <- fpf
  TPF <- tpf
  FPF <- FPF[ciIndx]
  TPF <- TPF[ciIndx]
  ciX <- binom.confint(x = FPF * K1, n = K1, method = "exact")
  ciY <- binom.confint(x = TPF * K2, n = K2, method = "exact")
  ciXUpper <- ciX$upper
  ciXLower <- ciX$lower
  ciYUpper <- ciY$upper
  ciYLower <- ciY$lower
  for (i in 1:length(FPF)){
    ciX <- data.frame(FPF = c(ciXUpper[i], ciXLower[i]), TPF = c(TPF[i], TPF[i]))
    ciY <- data.frame(FPF = c(FPF[i], FPF[i]), TPF = c(ciYUpper[i], ciYLower[i]))
    fitPlot <- fitPlot + geom_line(data = ciY, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = ciX, aes(x = FPF, y = TPF), color = "black")
    barRgt <- data.frame(FPF = c(ciXUpper[i], ciXUpper[i]), TPF = c(TPF[i] - 0.01, TPF[i] + 0.01))
    barLft <- data.frame(FPF = c(ciXLower[i], ciXLower[i]), TPF = c(TPF[i] - 0.01, TPF[i] + 0.01))
    barUp <- data.frame(FPF = c(FPF[i] - 0.01, FPF[i] + 0.01), TPF = c(ciYUpper[i], ciYUpper[i]))
    barBtm <- data.frame(FPF = c(FPF[i] - 0.01, FPF[i] + 0.01), TPF = c(ciYLower[i], ciYLower[i]))
    fitPlot <- fitPlot + geom_line(data = barRgt, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barLft, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barUp, aes(x = FPF, y = TPF), color = "black") + 
      geom_line(data = barBtm, aes(x = FPF, y = TPF), color = "black")
  }
  return(fitPlot)
}