# mainUnitNormalPdfCdf.R
rm(list = ls())
library(ggplot2)
x <- seq(-3,3,0.01)
pdf <- dnorm(x)
cdf <- pnorm(x)
df <- data.frame(
  z = c(x, x), all = c(pdf, cdf), 
  group = c(rep("pdf", length(pdf)), 
            rep("CDF", length(cdf)))
)

p <- ggplot(df, aes(x = z, y = all, color = group)) + 
  geom_line(data = df, size = 2) + xlab("z") + ylab("pdf/CDF") +
  geom_vline(xintercept = 1, linetype = 2, size = 2) + 
  scale_colour_manual(values=c("black","darkgrey")) + 
  theme(legend.title = element_blank(), legend.position = c(0.2, 0.95), 
        legend.key.size = unit(1.5, "lines"), legend.text=element_text(size=20, face = "bold"), 
        legend.direction = "horizontal") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))  +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 

print(p)

# qnorm(0.025)
# qnorm(1-0.025)
# pnorm(qnorm(0.025))
# qnorm(pnorm(-1.96))