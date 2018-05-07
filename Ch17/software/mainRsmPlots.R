# mainRsmPlots.R
rm(list = ls())
library(RJafroc)
library(ggplot2)

K2 <- 700;Lmax <- 4;Lk2 <- floor(runif(K2, 1, Lmax + 1))
nLesPerCase <- unique(Lk2);lesionDist <- array(dim = c(length(nLesPerCase), 2))
for (i in nLesPerCase) lesionDist[i, ] <- c(i, sum(Lk2 == i)/K2)

muArr <- c(2,3);lambda <- 1;nuArr <- c(0.15,0.25); L <- 1
for (i in 1:length(muArr)) {
  mu <- muArr[i]
  for (j in 1:length(nuArr)) {
    nu <- nuArr[j]
    ret1 <- PlotRsmOperatingCharacteristics(
      mu, lambda, nu, 
      type = "ALL", lesionDistribution = lesionDist, legendPosition  = "none"
    )
    pdfPlots <- ret1$PDFPlot + 
      scale_color_manual(values = c("black","darkgrey"), guide = FALSE) + 
      theme(axis.title.y = element_text(size = 25,face="bold"),
            axis.title.x = element_text(size = 30,face="bold"),
            legend.title = element_blank(), 
            legend.position = c(0.77,0.95), 
            legend.direction = "horizontal",
            legend.text = element_text(size = 20, face = "bold"),
            legend.key.size = unit(2, "lines")) +
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0)) 
    pdfPlots$layers[[1]]$aes_params$size <- 2 # line
    print(pdfPlots)
    cat("mu =", mu, ", lambda =", lambda,", nu =", nu,"\n")
  next
  }
}