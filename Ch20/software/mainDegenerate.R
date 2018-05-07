rm(list = ls()) # mainDegenerate.R
source("BMCurve.R")
source("CBMCurve.R")

plotOP <- data.frame(FPF = 0, TPF = 0.75)

a <- 0.6744898; b <- 0
plotCurve <- BMCurve(a, b)
figA <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = plotCurve) + 
  geom_point(data = plotOP) 
print(figA)

a <- 1.281552; b <- 0
plotCurve <- BMCurve(a, b)
figB <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = plotCurve) + 
  geom_point(data = plotOP)
print(figB)

a <- Inf; b <- 0
plotCurve <- BMCurve(a, b)
figC <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = plotCurve) + 
  geom_point(data = plotOP)
print(figC)

mu <- Inf; alpha <- 0.75
plotCurve <- CBMCurve(mu, alpha)
figD <- ggplot(mapping = aes(x = FPF, y = TPF)) + 
  geom_line(data = plotCurve) + 
  geom_point(data = plotOP)
print(figD)