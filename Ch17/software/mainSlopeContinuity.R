rm(list = ls()) # mainSlopeContinuity.R
library(RJafroc)
library(Rmpfr)

mu <- 0.5
lambda <- 0.1
nu <- 0.8
lambdaP <- lambda / mu
nuP <- 1 - exp(-mu * nu)
lesionDistr <- rbind(c(1, 0.2), c(2, 0.8))

for (myNegInf in (-3):(-15)) {
  ret <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, type = "ROC", 
    lesionDistribution = lesionDistr, 
    myNegInf = myNegInf) 
  
  cat("myNegInf = ", myNegInf, "\n")
  if (myNegInf == -3) print(ret$ROCPlot)
  
  i <- 1
  while (((ret$ROCPlot$data$TPF[i+1] - ret$ROCPlot$data$TPF[1]) == 0) ||
         ((ret$ROCPlot$data$FPF[i+1] - ret$ROCPlot$data$FPF[1]) == 0)){
    i <- i + 1
  }
  
  slopeContinuous <- (ret$ROCPlot$data$TPF[i+1]-ret$ROCPlot$data$TPF[1])/
    (ret$ROCPlot$data$FPF[i+1]-ret$ROCPlot$data$FPF[1])
  maxFPF <- (1 - exp(-lambdaP))
  maxTPF <- 1 - (0.2 * (1 - nuP)^1 * exp(-lambdaP) + 
    0.8 * (1 - nuP)^2 * exp(-lambdaP))
  slopeDashed <- (1 - maxTPF)/(1 - maxFPF)
  cat("slopeContinuous = ", slopeContinuous, ", slopeDashed = ", slopeDashed, "\n")
  
  zeta1 <- mpfr(-50, 2000) # set zeta1 = -50 with 2000 digit precision
  zeta2 <- zeta1 + 1e-12 # small increment
  delta_FPF <- (1 - exp(-lambdaP*pnorm(-zeta1))) - (1 - exp(-lambdaP*pnorm(-zeta2)))
  zeta <- zeta1
  TPF1 <- 0.2 * (1 - (1 - nuP * pnorm(mu-zeta))^1 * exp(-lambdaP*pnorm(-zeta))) + 
    0.8 * (1 - (1 - nuP * pnorm(mu-zeta))^2 * exp(-lambdaP*pnorm(-zeta)))
  zeta <- zeta2
  TPF2 <- 0.2 * (1 - (1 - nuP * pnorm(mu-zeta))^1 * exp(-lambdaP*pnorm(-zeta))) + 
    0.8 * (1 - (1 - nuP * pnorm(mu-zeta))^2 * exp(-lambdaP*pnorm(-zeta)))
  delta_TPF <- TPF1 - TPF2
  slopeNumeric <- delta_TPF/delta_FPF
  cat("slopeNumeric = ", as.numeric(slopeNumeric), "\n\n")
}
