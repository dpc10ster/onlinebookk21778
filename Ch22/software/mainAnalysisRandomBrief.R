# MainAnalysisRandomBrief.R
rm(list = ls()) 
library(RJafroc);library(ggplot2)
source("DiffFomAnal2007Hillis53.R");source("Wilcoxon.R")
# above line over-rides the RJAFROC version of DiffFomAnal2007Hillis53()
alpha <- 0.05
FOM <- "AUC" # allowed values are "PCL" "ALroc", "AUC" 
FPF <- 0.05
cat("FOM = ", FOM, "\n")
if (FOM == "PCL") cat("FPF = ", FPF, "\n")
cat("Random-reader random-case analysis")
cat("\nof Hupse Karssemeijer radiologist data:\n")
retNico <- DfReadLrocDataFile()
zjk1 <- retNico$NL[1,,,1]
zjk2Cl <- retNico$LLCl[1,,,1]
zjk2Il <- retNico$LLIl[1,,,1]
zjk2 <- pmax(zjk2Cl,zjk2Il)
K1 <- length(zjk1[1,])
K2 <- length(zjk2Cl[1,])
ret_nh <- DiffFomAnal2007Hillis53 (
  zjk1, 
  zjk2, 
  FOM, 
  FPF)  
cat(str(ret_nh))

