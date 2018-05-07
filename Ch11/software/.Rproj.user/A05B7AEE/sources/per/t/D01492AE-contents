rm(list = ls()) #mainNonCentralF.R  # freeze lines
library(ggplot2)
ndf <- 1;ddf <- 10;ncp <- c(0,2,5,10)
fCrit <- qf(0.95, ndf,ddf)
cat("critical value of x for rejecting NH is ", fCrit,"\n")
x <- seq(1, 20, 0.1)
for (i in 1:length(ncp))
{
  y <- df(x,ndf,ddf,ncp=ncp[i])
  cat("ndf = ", ndf, ", ddf = ", ddf, ", ncp = ", ncp[i], ", prob > fCrit = ", 
      1-pf(fCrit, ndf, ddf, ncp = ncp[i]), "\n")
  curveData <- data.frame(x = x, pdf = y)
  curvePlot <- ggplot(data = curveData, mapping = aes(x = x, y = pdf)) + geom_line(size = 1) +
    theme(axis.title.x = element_text(size = 30,face="bold"),
          axis.title.y = element_text(size = 30,face="bold"))
  print(curvePlot)
  stop("temp")
}