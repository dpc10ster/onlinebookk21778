# mainBestWorstObserver.R
rm(list = ls())
library(ggplot2)

frocDataPoints <- data.frame(NLF = c(0, 0, 1.45), LLF = c(1, 0, 0))
frocPlot <- ggplot(data = frocDataPoints, mapping = aes(x = NLF, y = LLF)) + 
  geom_line(size = 4) + 
  annotate("text", x = 0.1, y = 0.9, label = "Expert", size = 10) + 
  annotate("text", x = 0.75, y = 0.05, label = "Worst", size = 10) + 
  annotate("text", x = 0.75, y = 0.50, label = "Cannot define average FROC curve", size = 10) + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))
print(frocPlot)

afrocDataPoints <- data.frame(FPF = c(0, 0, 1), LLF = c(1, 0, 0))
afrocAvgPoints <- data.frame(FPF = c(0, 1), LLF = c(0.5, 0.5))
afrocDashed <- data.frame(FPF = c(0, 1, 1), LLF = c(1, 1, 0))
afrocPlot <- ggplot(mapping = aes(x = FPF, y = LLF)) + 
  geom_line(data = afrocDataPoints, size = 4) + 
  geom_line(data = afrocAvgPoints, size = 2) + 
  geom_line(data = afrocDashed, linetype = 2, size = 2) + 
  annotate("text", x = 0.1, y = 0.9, label = "Expert", size = 10) + 
  annotate("text", x = 0.5, y = 0.05, label = "Worst", size = 10) + 
  annotate("text", x = 0.5, y = 0.55, label = "Average", size = 10) + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))
print(afrocPlot)