rm(list = ls()) #mainNonCentralF.R 
library(ggplot2)
alpha  <- 0.05
ndf <- 1
ddf <- 100
ncpArr <- c(0,2,5,10)
fCrit <- qf(1-alpha, ndf,ddf)
rfCrit <- round(fCrit, 3)
cat("critical value for rejecting NH:", fCrit,"\n")
x <- seq(0, 20, 0.1)
for (i in 1:length(ncpArr))
{
  pdf <- df(x,ndf,ddf,ncp=ncpArr[i])
  cat("ndf = ", ndf, ", ddf = ", ddf, 
      "\nncp = ", ncpArr[i], "\nprob > fCrit = ", 
      1-pf(fCrit, ndf, ddf, ncp = ncpArr[i]), "\n")
  plotCurve <- data.frame(x = x, pdf = pdf)
  yMax <- 1 
  fText <- paste0("F[list(1-",as.character(alpha), 
                  ",", ndf, ",", ddf, ")] == ", rfCrit)
  p <- ggplot() + 
    geom_line(
      data = plotCurve, 
      mapping = aes(x = x, y = pdf), size = 1) + 
    xlab("z")
  print(p)
  next
}