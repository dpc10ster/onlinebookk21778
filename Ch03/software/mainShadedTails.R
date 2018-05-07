# MainShadedTails.R
rm(list = ls())
library(ggplot2)
mu <- 0;sigma <- 1
zeta <- -qnorm(0.025)
step <- 0.01

LL<- -3
UL <- mu + 3*sigma

x.values <- seq(zeta,UL,step)
cord.x <- c(zeta, x.values,UL) 
cord.y <- c(0,dnorm(x.values),0) 

z <- seq(LL, UL, by = step)
curveData <- data.frame(z = z, pdf = dnorm(z))
shadeData <- data.frame(z = cord.x, pdf = cord.y)
shadedTails <- ggplot(mapping = aes(x = z, y = pdf))  + 
  geom_polygon(data = shadeData, color = "grey", fill = "grey")

zeta <- qnorm(0.025)
x.values <- seq(LL, zeta,step)
cord.x <- c(LL, x.values,zeta) 
cord.y <- c(0,dnorm(x.values),0) 
shadeData <- data.frame(z = cord.x, pdf = cord.y)
shadedTails <- shadedTails + 
  geom_polygon(data = shadeData, color = "grey", fill = "grey") + 
  xlab("z")+ ylab("pdf" ) + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
shadedTails <- shadedTails + geom_line(data = curveData, color = "black", size = 2)
print(shadedTails)

