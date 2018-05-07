# mainRsmPdfRoc.R
rm(list = ls())
library(RJafroc)
library(ggplot2)

muArr <- c(0.001,1,2,3,4,5);lambda <- 1;nu <- 1; Lmax <- 1 
for (i in 1:length(muArr)) {
  mu <- muArr[i]
  ret1 <- PlotRsmOperatingCharacteristics(
    mu, lambda, nu, 
    type = "ALL", lesionDistribution = c(Lmax,1), legendPosition  = "none"
  )

  pdfPlots <- ret1$PDFPlot + 
    scale_color_manual(values = c("black","darkgrey"), guide = FALSE) + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold"),
          legend.title = element_blank(), 
          legend.position = c(0.25,0.95), 
          legend.direction = "horizontal", # use for D, E F
          legend.text = element_text(size = 20, face = "bold"),
          legend.key.size = unit(2, "lines")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0)) 
  pdfPlots$layers[[1]]$aes_params$size <- 2 # line
  print(pdfPlots)
  #next
  
  ROCPlot <- ret1$ROCPlot + scale_color_manual(values = "black")
  fpfMax <- max(ret1$ROCPlot$data$FPF)
  tpfMax <- max(ret1$ROCPlot$data$TPF)
  if (fpfMax < 0.99){
    fpfCross <- (fpfMax + tpfMax) / 2
    tpfCross <- fpfCross
    endPoint <- data.frame(FPF = fpfMax, TPF = tpfMax)
    ds <- data.frame(FPF = c(fpfMax, fpfCross), TPF = c(tpfMax, tpfCross))
    diagonal <- data.frame(FPF = c(0, 1), TPF = c(0, 1))
    dsText <- data.frame(FPF = (fpfMax + fpfCross)/2 + 0.05, TPF = (tpfMax + tpfCross)/2)
    ROCPlot <- ROCPlot +
      geom_point(data = endPoint, 
                 mapping = aes(x = FPF, y = TPF), 
                 shape = 15, size = 7) +
      geom_line(data = ds, 
                mapping = aes(x = FPF, y = TPF), 
                linetype = 2, size = 2) +
      geom_line(data = diagonal, 
                mapping = aes(x = FPF, y = TPF), 
                linetype = 2, size = 2) +
      geom_text(data = dsText, 
                mapping = aes(x = FPF, y = TPF), 
                label = "d[s]", parse = TRUE, size = 10)
  }
  ROCPlot <- ROCPlot + 
    theme(axis.title.y = element_text(size = 25,face="bold"),
          axis.title.x = element_text(size = 30,face="bold")) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(expand = c(0, 0))
  ROCPlot$layers[[1]]$aes_params$size <- 2 # line
  ROCPlot$layers[[2]]$aes_params$size <- 2 # line
  print(ROCPlot) 
  cat("mu = ", mu,
      "\nlambda = ", lambda,
      "\nnu = ", nu, 
      "\nAUC = ", ret1$aucROC,
      "\nfpfMax = ", fpfMax,
      "\ntpfMax = ", tpfMax,"\n")
}