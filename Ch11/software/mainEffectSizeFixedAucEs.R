rm(list = ls()) #mainEffectSizeFixedAucEs.R
source("AzDeePrimeTransformations.R")
require(ggplot2)

esAuc = 0.05
nBins <- 100
AzNh  <-  seq(0.55,1,length.out = nBins)
esDp <-  effectSizeDeePrime(esAuc, AzNh)
AzNh <-  AzNh[1:length(esDp)]
esDpRatio  <-  esDp / azToDeePrime(AzNh)
ConstAucPlot  <-  rep(esAuc, length.out = length(esDpRatio))
df <- data.frame(Az = c(AzNh, AzNh),DpMultipler = c(ConstAucPlot, esDpRatio), 
                 truth = c(rep('ES: dp multiple', length(AzNh)), rep('Const. Az ES', length(AzNh))))

myPlot <- ggplot(df, aes(x = Az, y = DpMultipler, color = truth)) + 
  geom_line(size = 1) +
  scale_colour_manual(values=c("black","darkgrey")) + 
  theme(legend.title = element_blank(), legend.position = c(0.52, 0.85)) +
  theme(axis.title.x = element_text(size = 30,face="bold"),
        axis.title.y = element_text(size = 30,face="bold")) +
  theme(legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(2.5, "lines"))


print(myPlot)

