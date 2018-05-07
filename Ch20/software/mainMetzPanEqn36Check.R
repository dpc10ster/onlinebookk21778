# mainMetzPanEqn36Check.R
rm( list = ls()) 
require(ggplot2)
require(caTools)
require(mvtnorm)
source("proprocFunctions.R")
source("TonyData.R") # this contains the saved values from the PROPROC run

npts <-  10000
for (i in 1:2) {
  for (j in 1:5) {
    C  <-  c1[i,j]
    da  <-  d_a1[i,j]
    ret <- GetLimits(da,C)
    LL <- ret$LL;UL <- ret$UL
    vc  <-  seq (LL, UL, length.out = npts)
    TPF  <-  TruePositiveFraction (vc, da, C)
    FPF <- FalsePositiveFraction (vc, da, C)
    FPF <- rev(FPF);TPF <- rev(TPF)
    df2 <- data.frame(FPF = FPF, TPF = TPF)
    plotRoc <- ggplot(df2, aes(x = FPF, y = TPF)) + geom_line()
#    print(plotRoc) not to be overwhelmed with plots

    # do integral numerically
    numAuc <- trapz(FPF, TPF)
    
    # Implement Eqn. 36 from Metz-Pan paper 
    rho <- -(1-C^2)/(1+C^2);sigma <- rbind(c(1, rho), c(rho, 1))
    lower <- rep(-Inf,2);upper <- c(-da/sqrt(2),0)
    A_prop <- pnorm(da/sqrt(2)) + 2 * pmvnorm(lower, upper, sigma = sigma)
    A_prop <-  as.numeric(A_prop)
    
    cat("i = ", i,"j = ", j,"C = ", C, ", da = ", da, "NumericalAUC = ", numAuc, ", Eqn. 36 = ", A_prop,"\n")
  }
}