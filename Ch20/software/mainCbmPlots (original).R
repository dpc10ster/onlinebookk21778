rm(list = ls()) #mainCbmRoc.R
source("cbmFunctions.R")

FPF <- seq(0.0, 1, 0.001)
alpha <- 0.8;mu <- 3
TPF <- CbmRocY(FPF, mu, alpha)

rocPlot <- data.frame(FPF = FPF, TPF = TPF)
plotRoc <- ggplot(rocPlot, aes(x = FPF, y = TPF)) + geom_line()
print(plotRoc)

z <- seq(-3, 5, by = 0.01) # may need to adjust limits to view detail of slope plot

slope <- ((1-alpha)*dnorm(-z) + alpha*dnorm(mu-z))/dnorm(-z) # same as likelihood ratio

slopePlot <- data.frame(z = z, slope = slope)
plotSlope <- ggplot(slopePlot, aes(x = z, y = slope)) + geom_line()
print(plotSlope)

z1 <- seq(-3, 3, by = 0.01)
z2 <- seq(-3, mu + 3, by = 0.01)

Pdf1 <- dnorm(z1)
Pdf2 <- (1 - alpha) * dnorm(z2) + alpha * dnorm(z2, mu)

df <- data.frame(
  z = c(z1, z2), pdf = c(Pdf1, Pdf2), 
  truth = c(rep('non-diseased', length(Pdf1)), 
            rep('diseased', length(Pdf2)))
)

cbmPdfs <- ggplot(df, aes(x = z, y = pdf, color = truth)) + 
  geom_line() + 
  scale_colour_manual(values=c("red","green")) + 
  theme(legend.title = element_blank(), legend.position = c(0.9, 0.9))

print(cbmPdfs)
cat("mu = ", mu, ", alpha = ", alpha, "\n")

