# mainRocPdfs.R
rm( list = ls())
require(ggplot2)

mu <- 1.5
sigma <- 1.5

z1 <- seq(-3, 3, by = 0.01)
z2 <- seq(-3, 7, by = 0.01)

pdf1 <- dnorm(z1)
pdf2 <- dnorm(z2, mu, sd = sigma)

df <- data.frame(
  z = c(z1, z2), 
  pdf = c(pdf1, pdf2), 
  truth = c(rep('non-diseased', length(pdf1)), 
  rep('diseased', length(pdf2))))

rocPdfs <- ggplot(
  df, 
  aes(x = z, y = pdf, color = truth)) + 
  geom_line() + 
  scale_colour_manual(values=c("red","green")) + 
  theme(legend.title = element_blank(), legend.position = c(0.9, 0.9))

print(rocPdfs)
