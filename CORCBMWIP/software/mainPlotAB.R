rm( list = ls())
library(ggplot2)

mu <- 2
sigma <- 1.8
x <- seq(-3, 6, by = 0.01)
y1 <- dnorm(x)
y2 <- dnorm(x, mu, sigma)

y <- data.frame(x = c(x, x), y = c(y1, y2), curves = c(rep("y1", length(y1)), rep("y2", length(y2))))

p1 <- ggplot() + geom_line(data = y, mapping = aes(x = x, y = y, linetype = curves)) + 
  scale_linetype_manual(values = c("dashed", "solid")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black"), legend.position="none",
        axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) + 
  geom_segment(aes(x = 1, y = 0.405, xend = 0, yend = 0.405), arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 1, y = 0.405, xend = mu, yend = 0.405), arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 1, y = 0.39, label = "a", size = 10) + 
  annotate("text", x = 1.2, y = 0.35, label = "Noise", size = 10) + 
  geom_segment(aes(x = 0, y = 0.205, xend = 1.15, yend = 0.205), arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 0, y = 0.205, xend = -1.15, yend = 0.205), arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 0, y = 0.19, label = "b", size = 10) + 
  geom_segment(aes(x = 2, y = 0.105, xend = -0.2, yend = 0.105), arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 2, y = 0.105, xend = 4.2, yend = 0.105), arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 2, y = 0.09, label = "1", size = 10) + 
  annotate("text", x = 4.2, y = 0.18, label = "Signal", size = 10) 
print(p1)

p2 <- ggplot() + geom_line(data = y, mapping = aes(x = x, y = y, linetype = curves)) +
  scale_linetype_manual(values = c("dashed", "solid")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black"), legend.position="none",
        axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) + 
  geom_segment(aes(x = 1, y = 0.405, xend = 0, yend = 0.405), arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 1, y = 0.405, xend = mu, yend = 0.405), arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 1, y = 0.39, label = "mu == a/b", size = 10, parse = TRUE) + 
  annotate("text", x = 1.2, y = 0.35, label = "Noise", size = 10) + 
  geom_segment(aes(x = 0, y = 0.205, xend = 1.15, yend = 0.205), arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 0, y = 0.205, xend = -1.15, yend = 0.205), arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 0, y = 0.19, label = "1", size = 10) + 
  geom_segment(aes(x = 2, y = 0.105, xend = -0.2, yend = 0.105), arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 2, y = 0.105, xend = 4.2, yend = 0.105), arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 2, y = 0.09, label = "sigma == 1/b", size = 10, parse = TRUE) + 
  annotate("text", x = 4.2, y = 0.18, label = "Signal", size = 10) 
print(p2)
