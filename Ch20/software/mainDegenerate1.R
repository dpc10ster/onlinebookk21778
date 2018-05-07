rm(list = ls()) # mainDegenerate.R
source("BMCurve.R")
source("CBMCurve.R")

plotOP <- data.frame(FPF = 0, TPF = 0.75)

a <- 0.6744898; b <- 0
plotCurve <- BMCurve(a, b)
figA <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = plotCurve, size = 4) + 
  geom_point(data = plotOP, size = 10)  + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
print(figA)

a <- 1.281552; b <- 0
plotCurve <- BMCurve(a, b)
figB <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = plotCurve, size = 4) + 
  geom_point(data = plotOP, size = 10)  + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
print(figB)

a <- Inf; b <- 0
plotCurve <- BMCurve(a, b)
figC <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = plotCurve, size = 4) + 
  geom_point(data = plotOP, size = 10)  + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
print(figC)

mu <- Inf; alpha <- 0.75
plotCurve <- CBMCurve(mu, alpha)
figD <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = plotCurve, size = 4) + 
  geom_point(data = plotOP, size = 10)  + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
print(figD)