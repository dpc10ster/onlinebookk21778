rm( list = ls()) # mainRocSlopeTonyData.R
source("ProprocFunctions.R")
require(ggplot2)
require(caTools)
require(mvtnorm)
source("TonyData.R")

npts <-  10000
for (i in 1:2) {
  for (j in 1:5) {
    C  <-  c1[i,j]
    da  <-  d_a1[i,j]
    ret <- Transform2ab(da, C)
    a <- ret$a;b <- ret$b
    z <- seq(-3, 5, by = 0.01) # may need to adjust limits to view detail of slope plot
    ret <- GetLimits(da,C)
    LL <- ret$LL;UL <- ret$UL
    vc  <-  seq (LL, UL, length.out = npts)
    TPF  <-  TruePositiveFraction (vc, da, C)
    FPF <- FalsePositiveFraction (vc, da, C)
    FPF <- rev(FPF);TPF <- rev(TPF)
    df2 <- data.frame(FPF = FPF, TPF = TPF)
    plotRoc <- ggplot(df2, aes(x = FPF, y = TPF)) + geom_line()
    print(plotRoc)

    slope <-b*dnorm(a-b*z)/dnorm(-z) # same as likelihood ratio
    
    slopePlot <- data.frame(z = z, slope = slope)
    plotSlope <- ggplot(slopePlot, aes(x = z, y = slope)) + geom_line()
    print(plotSlope)
    next
  }
}
