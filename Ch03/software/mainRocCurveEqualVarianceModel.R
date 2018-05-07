# MainRocCurveEqualVarianceModel.R
rm(list = ls())
source("drawROC.R")
library(ggplot2)
mu <- 0
zeta <- seq(-5, mu + 5, 0.05)
FPF <- pnorm(-zeta)
rocPlot <- ggplot(mapping = aes(x = FPF, y = TPF))
for (mu in 0:3){
  TPF <- pnorm(mu-zeta)
  curveData <- data.frame(FPF = FPF, TPF = TPF)
  rocPlot <- rocPlot + 
    geom_line(data = curveData, size = 2) + 
    xlab("FPF")+ ylab("TPF" ) + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"))  +
  annotate("text", x = pnorm(-mu/2) + 0.07, y = pnorm(mu/2), 
             label = paste0("mu == ", mu), parse = TRUE, size = 8)
  next
}
rocPlot <- rocPlot +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))     
  
rocPlot <- rocPlot + 
  geom_abline(slope = -1, intercept = 1, linetype = 3, size = 2)
print(rocPlot)