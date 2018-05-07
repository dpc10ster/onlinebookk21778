rm(list = ls()) #mainEffectSizeFixedDpMultiple.R
source("AzDeePrimeTransformations.R")
require(ggplot2)

nBins <- 100
dpMultiple = 0.2
azNh  <-  seq(0.55,0.99,length.out = nBins)
dpAH <- azToDeePrime(azNh) * (1 + dpMultiple)
azAh <-  deePrimeToAz(dpAH)
esAz <- azAh - azNh
ConstDpPlot  <-  rep(dpMultiple, length.out = nBins)
df <- data.frame(Az = c(azNh, azNh),esAuc = c(esAz, ConstDpPlot), 
                 truth = c(rep('Az ES for fixed dp ES', length(azNh)), rep('ES: dp multiplier', length(azNh))))

myPlot <- ggplot(df, aes(x = Az, y = esAuc, color = truth)) + 
  geom_line(size = 1) +
  scale_colour_manual(values=c("black","darkgrey")) + 
  theme(legend.title = element_blank(), legend.position = c(0.52, 0.7)) +
  theme(axis.title.x = element_text(size = 30,face="bold"),
        axis.title.y = element_text(size = 30,face="bold")) +
  theme(legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(2.5, "lines"))

print(myPlot)