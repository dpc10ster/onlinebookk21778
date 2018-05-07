rm(list = ls()) #mainFDist.R
require(ggplot2)
alpha  <- 0.05;ndf <- 1;ddf <- 200
fCrit <- qf(1 - alpha, ndf, ddf)
rfCrit <- round(fCrit, 3)
cat("alpha = ", alpha, "ndf = ", ndf, "ddf = ", 
    ddf, "fCrit = ", fCrit, "\n")

x <- seq(0, 20, by = 0.01)
pdf <- df(x, ndf, ddf)
plotCurve <- data.frame(x = x, pdf = pdf)
plotArea <- data.frame(
  x = c(0, seq(0, fCrit, by = 0.01), fCrit), 
  y = c(0, df(seq(0, fCrit, by = 0.01), 
              ndf, ddf), 0))
yMax <- 1 
fText <- paste0(
  "F[list(1-",
  as.character(alpha), 
  ",", ndf, ",", ddf, ")] == ", rfCrit)
p <- ggplot() + 
  geom_line(
    data = plotCurve, 
    mapping = aes(x = x, y = pdf), size = 2) + 
  geom_polygon(
    data = plotArea, 
    mapping = aes(x = x, y = y)) +
  geom_segment(
    aes(x = fCrit + 3, 
        y = yMax / 4, 
        xend = fCrit, 
        yend = df(fCrit, ndf, ddf)), 
    arrow = arrow(
      length = unit(0.5, "cm")), size = 2) +
  annotate("text", x = rfCrit + 4, y = yMax / 4 + 0.04, 
           label = fText, size = 7, parse = TRUE) + 
  coord_cartesian(ylim = c(0, yMax))

print(p)