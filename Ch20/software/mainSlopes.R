# mainSlopes.R
rm( list = ls())
require(ggplot2)

# values used in Fig. 20.4 (A-F) of book
a <- 0.7;b <- 0.5 
zArray <- list(seq(-5, 8, by = 0.01), 
               seq(-3, 5, by = 0.01), 
               seq(-3, 3, by = 0.01), 
               seq(-2, 2, by = 0.01), 
               seq(-1, 1, by = 0.01), 
               seq(-1, 0, by = 0.01))

for (i in 1:length(zArray)) 
{
  z <- zArray[[i]]
  slope <-b*dnorm(a-b*z)/dnorm(-z) # same as likelihood ratio
  slopePlot <- data.frame(z = z, slope = slope)
  plotSlope <- ggplot(slopePlot, aes(x = z, y = slope)) + 
    geom_line(size = 2)  + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"))  +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  print(plotSlope)
  cat("a = ", a, ", b = ", b, "\n")
}