rm(list = ls()) # mainAnalysisFixed.R; i.e., case regarded as a fixed effect
library("RJafroc") 
library(ggplot2)
#require("caTools") #needed for trapezoidal area

alpha <- 0.05
cat("Hupse-Karssemeijer analysis\n: random readers fixed cases\n")

FOM <- "Wilcoxon" # allowed values are "PCL" "ALROC", "Wilcoxon" 
FPFArr <- c(0.05, 0.2, 0.5, 1) # at which to evaluate PCL

zjk1 <- datasetCadLroc$NL[1,,,1];zjk2Cl <- datasetCadLroc$LLCl[1,,,1];zjk2Il <- datasetCadLroc$LLIl[1,,,1];zjk2 <- pmax(zjk2Cl,zjk2Il);J <- dim(zjk1)[1] - 1
for (i in 1:length(FPFArr)) {
  FPF <- FPFArr[i]
  cat("\nFOM = ", FOM, "\n");if (FOM != "Wilcoxon") cat("FPF = ", FPF, "\n")
  
  thetajc <- array (dim = (J+1))
  if (FOM == "Wilcoxon") {
    for (j in 1:(J+1)) thetajc[j] <- Wilcoxon(zjk1[j,],zjk2[j,]) 
  } else if (FOM == "PCL") {
    for (j in 1:(J+1)) thetajc[j] <- (LrocFoms(zjk1[j,], zjk2Cl[j,], FPF))$PCL 
  } else if (FOM == "ALROC") {
    for (j in 1:(J+1)) thetajc[j] <- (LrocFoms(zjk1[j,], zjk2Cl[j,], FPF))$ALroc
  } else stop("wrong FOM value")  
  
  fom_diff <- thetajc[-1] - thetajc[1]
  ret <- t.test(fom_diff)
  
  cat("FomCad = ", thetajc[1],"\n")
  cat("AVG radiologist performance = ", mean(thetajc[-1]),
      "\n95%CI = ", as.numeric(ret$conf.int)+thetajc[1],"\n")
  cat("AVG diff. performance = ", mean(thetajc[-1]-thetajc[1]),
      "\n95%CI = ", as.numeric(ret$conf.int),"\n")
  cat("t-statistic = ", as.numeric(ret$statistic),"\n")
  cat("df = ", as.numeric(ret$parameter),"\n")
  cat("p-value = ", as.numeric(ret$p.value),"\n\n")
  if ((i == 1) && (FOM != "Wilcoxon")) {
    plots1 <- LrocPlots (zjk1, zjk2Cl, seq(1,9))$lrocPlot # last argument is the readers to display, all 9 readers in this case  
    plots1 <- plots1 + 
      theme(legend.title = element_blank(), legend.position = c(0.9, 0.05), 
            legend.key.size = unit(1.5, "lines"), legend.text=element_text(size=20, face = "bold")) + 
      theme(axis.title.y = element_text(size = 25,face="bold"),
            axis.title.x = element_text(size = 30,face="bold"))  +
      scale_x_continuous(limits = c(0,1)) + 
      scale_y_continuous(limits = c(0,1)) 
    print(plots1)
  }
}

