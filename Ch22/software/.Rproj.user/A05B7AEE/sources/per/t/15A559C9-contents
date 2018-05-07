# MainAnalysisRandom.R # local version of functions; see also CheckInterpolation.xlsx
rm(list = ls()) 
library(RJafroc);library(ggplot2);library("caTools") #needed for trapezoidal area
source("DiffFomAnal2007Hillis53.R");source("Wilcoxon.R")
source("LrocFoms1.R");source("LrocOperatingPointsFromRatings.R")

cat("Random-reader random-case analysis")
cat("\nof Hupse Karssemeijer radiologist data:\n")
FOM <- "ALROC" # allowed values are "PCL" "ALROC", "Wilcoxon" 
FPFArr <- c(0.05, 0.2, 0.5, 1)
retNico <- DfReadLrocDataFile()
zjk1 <- retNico$NL[1,,,1]
zjk2Cl <- retNico$LLCl[1,,,1]
zjk2Il <- retNico$LLIl[1,,,1]
zjk2 <- pmax(zjk2Cl,zjk2Il)
for (i in 1:length(FPFArr)) {
  FPF <- FPFArr[i]
  cat("\n\nFOM = ", FOM, "\n")
  if (FOM != "Wilcoxon") cat("FPF = ", FPF, "\n")
  if (FOM == "Wilcoxon") {
    ret <- DiffFomAnal2007Hillis53 (
      zjk1, 
      zjk2, 
      FOM)  
  } else if (FOM == "PCL") {
    ret <- DiffFomAnal2007Hillis53 (
      zjk1, 
      zjk2Cl, 
      FOM, 
      FPF)  
  } else if (FOM == "ALROC") {
    ret <- DiffFomAnal2007Hillis53 (
      zjk1, 
      zjk2Cl, 
      FOM, 
      FPF)  
  } else stop("wrong FOM value")
  print(ret)
  thetajc <- ret$thetajc
  cat("FOM Cad = ", ret$thetajc[1], "\n")
  cat("FOM Rad = ", mean(ret$thetajc[-1]), "\n")
  psijc <- thetajc[-1] - thetajc[1]
  avgRad <- mean(ret$thetajc[-1])
  sdRad  <- sd(psijc)
  CIRad <- avgRad + ret$CI - ret$PsiMean
  cat("CIRad = ", CIRad, "\n")
  if (FOM == "Wilcoxon") break
}