#MainShadedPlots.R
rm(list = ls())
library(ggplot2)
mu <- 3;sigma <- 1;zeta <- 1;options(digits=3)
step <- 0.05

lowerLimit<- -1 # lower limit
upperLimit <- mu + 3*sigma # upper limit

z <- seq(lowerLimit, upperLimit, by = step)
pdfs <- dnorm(z)
seqNor <- seq(zeta,upperLimit,step)
cord.x <- c(zeta, seqNor,upperLimit) 
# need two y-coords at each end point of range; 
# one at zero and one at value of function
cord.y <- c(0,dnorm(seqNor),0) 
curveData <- data.frame(z = z, pdfs = pdfs)
shadeData <- data.frame(z = cord.x, pdfs = cord.y)
shadedPlots <- ggplot(mapping = aes(x = z, y = pdfs)) + 
  geom_line(data = curveData, color = "blue", size = 2) + 
  geom_polygon(data = shadeData, color = "blue", fill = "blue")

crossing <- uniroot(function(x) dnorm(x) - dnorm(x,mu,sigma), 
                    lower = 0, upper = 3)$root
crossing <- max(c(zeta, crossing))
seqAbn <- seq(crossing,upperLimit,step)
cord.x <- c(seqAbn, rev(seqAbn))
# reason for reverse 
# we want to explicitly define the polygon
# we dont want R to close it 

cord.y <- c()
for (i in seq(1,length(cord.x)/2)) {
  cord.y <- c(cord.y,dnorm(cord.x[i],mu, sigma))
}
for (i in seq(1,length(cord.x)/2)) {
  cord.y <- c(cord.y,dnorm(cord.x[length(cord.x)/2+i]))
}
pdfs <- dnorm(z, mu, sigma)
curveData <- data.frame(z = z, pdfs = pdfs)
shadeData <- data.frame(z = cord.x, pdfs = cord.y)
shadedPlots <- shadedPlots + geom_line(data = curveData, color = "red", size = 2) + 
  geom_polygon(data = shadeData, color = "red", fill = "red")
seqAbn <- seq(zeta,upperLimit,step)
for (i in seqAbn) {
  # define xs and ys of two points, separated only along y-axis
  vlineData <- data.frame(x1 = i, x2 = i, y1 = 0, y2 = dnorm(i, mu, sigma))
  # draw vertical line between them
  shadedPlots <- shadedPlots + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), data = vlineData, color = "red")
}
shadedPlots <- shadedPlots + xlab("z")+ ylab("pdfs" ) + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
print(shadedPlots)