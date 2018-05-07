rm(list = ls()) #mainNonCentralF.R  # freeze lines
library(ggplot2)
alpha  <- 0.05;ndf <- 2;ddf <- 100;ncp <- c(0,2,5,10)
fCrit <- qf(1-alpha, ndf,ddf)
rfCrit <- round(fCrit, 3)
cat("critical value of x for rejecting NH is ", fCrit,"\n")
x <- seq(0, 20, 0.01)
for (i in 1:length(ncp))
{
  pdf <- df(x,ndf,ddf,ncp=ncp[i])
  cat("ndf = ", ndf, ", ddf = ", ddf, ", ncp = ", ncp[i], ", prob > fCrit = ", 
      1-pf(fCrit, ndf, ddf, ncp = ncp[i]), "\n")
  plotCurve <- data.frame(x = x, pdf = pdf)
  plotArea <- data.frame(x = c(0, seq(0, fCrit, by = 0.01), fCrit), 
                         y = c(0, df(seq(0, fCrit, by = 0.01), ndf, ddf,ncp=ncp[i]), 0))
  yMax <- 1 
  fText <- paste0("F[list(1-",as.character(alpha), ",", ndf, ",", ddf, ")] == ", rfCrit)
  p <- ggplot() + geom_line(data = plotCurve, mapping = aes(x = x, y = pdf), size = 1) + 
    xlab("z") +
    theme(axis.title.x = element_text(hjust = 0.8, size = 30,face="bold"),
          axis.title.y = element_text(size = 30,face="bold")) +
    # geom_polygon(data = plotArea, mapping = aes(x = x, y = y)) +
    # geom_segment(aes(x = fCrit + 3, y = yMax / 4, xend = fCrit, yend = df(fCrit, ndf, ddf,ncp=ncp[i])), 
    #              arrow = arrow(length = unit(0.3, "cm"))) +
    # annotate("text", x = rfCrit + 4, y = yMax / 4 + 0.04, 
    #          label = fText, size = 7, parse = TRUE) + 
    # coord_cartesian(ylim = c(0, yMax)) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  
  print(p)
  next
}