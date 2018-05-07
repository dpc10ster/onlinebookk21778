# this version uses the RJafroc version of DiffFomAnal2007Hillis53()
# MainAnalysisRandom1.R
rm(list = ls()) 
library(RJafroc)
library(ggplot2)
library("caTools") #needed for trapezoidal area
cat("Random-reader random-case analysis")
cat("\nof Hupse Karssemeijer radiologist data:\n")
FOM <- "ALROC" # allowed values are "PCL" "ALROC", "Wilcoxon" 
FPFArr <- c(0.05, 0.2, 0.5, 1)
retNico <- DfReadLrocDataFile()
for (i in 1:length(FPFArr)) {
  FPF <- FPFArr[i]
  cat("FOM = ", FOM, "\n")
  if (FOM == "PCL") cat("FPF = ", FPF, "\n")
  if (FOM == "Wilcoxon") {
    ret_nh2 <- DiffFomAnal2007Hillis53 (
      datasetCadLroc,
      FOM = FOM)
  } else if (FOM == "PCL") {
    ret_nh2 <- DiffFomAnal2007Hillis53 (
      retNico,
      FOM = FOM,
      FPFValue = FPF)
  } 
  else if (FOM == "ALROC") {
    ret_nh2 <- DiffFomAnal2007Hillis53 (
      retNico,
      FOM = FOM,
      FPFValue = FPF)
  } else stop("wrong FOM value")
  print(ret_nh2)
  if (FOM == "Wilcoxon") break
}