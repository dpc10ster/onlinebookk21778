#MainPropertiesUnitNormal.R
rm(list = ls());library(ggplot2)
cat("2.5 percentile of distribution =", qnorm(0.025), "\n")
cat("97.5 percentile of distribution =", qnorm(1-0.025), "\n")
cat("P(X<0) =", pnorm(0),"\n")
cat("P(X<-1.96) =", pnorm(-1.96),"\n")
cat("P(X<-Infinity) =", pnorm(-Inf),"\n")
cat("P(X<Infinity) =", pnorm(Inf),"\n")
# plot the CDF
x <- seq(-3, 3, by = 0.01)
curveData <- data.frame(x = x, CDF = pnorm(x))
cdfPlot <- ggplot(data = curveData, mapping = aes(x = x, y = CDF)) + geom_line()
print(cdfPlot)