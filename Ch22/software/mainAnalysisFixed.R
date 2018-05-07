# mainAnalysisFixed.R
# case regarded as a fixed effect
rm(list = ls()) 
library("RJafroc") 
library(ggplot2)
require("caTools") #needed for trapezoidal area
source("Wilcoxon.R");source("LrocFoms1.R");source("LrocOperatingPointsFromRatings.R")
alpha <- 0.05
cat("Hupse-Karssemeijer analysis\n: random readers fixed cases\n")

FOM <- "PCL" # allowed values are "PCL" "ALROC", "Wilcoxon" 
FPFArr <- c(0.2, 0.5, 1) # at which to evaluate PCL #c(0.05, 0.2, 0.5, 1)

retNico <- DfReadLrocDataFile()
zjk1 <- retNico$NL[1,,,1]
zjk2Cl <- retNico$LLCl[1,,,1]
zjk2Il <- retNico$LLIl[1,,,1]
zjk2 <- pmax(zjk2Cl,zjk2Il)
J <- dim(zjk1)[1] - 1
for (i in 1:length(FPFArr)) {
  FPF <- FPFArr[i]
  cat("\n\nFOM = ", FOM, "\n")
  if (FOM != "Wilcoxon") cat("FPF = ", FPF, "\n")
  
  thetajc <- array (dim = (J+1))
  if (FOM == "Wilcoxon") {
    for (j in 1:(J+1)) 
      thetajc[j] <- Wilcoxon(zjk1[j,],zjk2[j,]) # note subtle diff. in 2nd argument
  } else if (FOM == "PCL") {
    for (j in 1:(J+1)) 
      thetajc[j] <- (LrocFoms(zjk1[j,], zjk2Cl[j,], FPF))$PCL # note subtle diff. in 2nd argument 
  } else if (FOM == "ALROC") {
    for (j in 1:(J+1)) 
      thetajc[j] <- (LrocFoms(zjk1[j,], zjk2Cl[j,], FPF))$ALroc # extracting ALROC instead of PCL in line 31
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
  if ((i == 1) && (FOM == "PCL")) {
    # last argument is the readers to display, all 9 readers in this case
    plots1 <- LrocPlots (zjk1, zjk2Cl, seq(1,9))$lrocPlot   
    print(plots1)
  }
  if (FOM == "Wilcoxon") break
}

