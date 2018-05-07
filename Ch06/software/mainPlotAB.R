rm( list = ls())# mainPlotAB.R
library(ggplot2) 

mu <- 2
sigma <- 1.8;a <- mu/sigma;b <- 1/sigma
z <- seq(-3, 6, by = 0.01)
y1 <- dnorm(z)
y2 <- dnorm(z, mu, sigma)

y <- data.frame(z = c(z, z), pdfs = c(y1, y2), curves =
                  c(rep("y1", length(y1)), rep("y2", length(y2))))

p1 <- ggplot() + geom_line(data = y, mapping = aes(x = z, y = pdfs, linetype = curves), size = 2) +
  scale_linetype_manual(values = c("dashed", "solid")) + theme_bw() +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black"), legend.position="none") +
  geom_segment(aes(x = 1, y = 0.405, xend = 0, yend = 0.405),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 1, y = 0.405, xend = mu, yend = 0.405),
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = 1, y = 0.39, label = "a", size = 10, fontface = "bold") +
  annotate("text", x = 1.2, y = 0.35, label = "Noise", size = 10, fontface = "bold") +
  geom_segment(aes(x = 0, y = 0.205, xend = 1.15, yend = 0.205),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 0, y = 0.205, xend = -1.15, yend = 0.205),
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = 0, y = 0.19, label = "b", size = 10, fontface = "bold") +
  geom_segment(aes(x = 2, y = 0.105, xend = -0.2, yend = 0.105),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 2, y = 0.105, xend = 4.2, yend = 0.105),
               arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = 2, y = 0.09, label = "1", size = 10, fontface = "bold") +
  annotate("text", x = 4.2, y = 0.18, label = "Signal", size = 10, fontface = "bold")
print(p1)

z <- seq(-3*sigma, 6*sigma, by = 0.01)
y1 <- dnorm(z, sd = sigma)
y2 <- dnorm(z, mean = mu*sigma, sd = sigma^2)

y <- data.frame(z = c(z, z), pdfs = c(y1, y2), curves = 
                  c(rep("y1", length(y1)), rep("y2", length(y2))))
p2 <- ggplot() + geom_line(data = y, mapping = aes(x = z, y = pdfs, linetype = curves), size = 2) +
  scale_linetype_manual(values = c("dashed", "solid")) + theme_bw() +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black"), legend.position="none") + 
  geom_segment(aes(x = 1 * sigma, y = 0.405 / sigma, xend = 0, yend = 0.405 / sigma), 
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 1 * sigma, y = 0.405 / sigma, xend = mu * sigma, yend = 0.405 / sigma), 
               arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 1 * sigma, y = 0.39 / sigma, label = "mu == a/b", 
           size = 10, fontface = "bold", parse = TRUE) + 
  annotate("text", x = 1.2 * sigma, y = 0.35 / sigma, label = "Noise", 
           size = 10, fontface = "bold", parse = TRUE) + 
  geom_segment(aes(x = 0, y = 0.205 / sigma, xend = 1.15 * sigma, yend = 0.205 / sigma), 
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 0, y = 0.205 / sigma, xend = -1.15 * sigma, yend = 0.205 / sigma), 
               arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 0, y = 0.19 / sigma, label = "1", size = 10, fontface = "bold") + 
  geom_segment(aes(x = 2 * sigma, y = 0.105 / sigma, xend = -0.2 * sigma, yend = 0.105 / sigma), 
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_segment(aes(x = 2 * sigma, y = 0.105 / sigma, xend = 4.2 * sigma, yend = 0.105 / sigma), 
               arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text", x = 2 * sigma + 0.1, y = 0.09 / sigma, label = "sigma == 1/b", 
           size = 10, fontface = "bold", parse = TRUE) + 
  annotate("text", x = 4.2 * sigma, y = 0.18 / sigma, label = "Signal", 
           size = 10, fontface = "bold", parse = TRUE)
  print(p2)
