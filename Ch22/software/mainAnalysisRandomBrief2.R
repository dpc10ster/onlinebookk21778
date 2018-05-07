# MainAnalysisRandomBrief2.R
rm(list = ls()) 
library(RJafroc)
library(ggplot2)
# this version uses the RJafroc version of DiffFomAnal2007Hillis53()
alpha <- 0.05
FOM <- "Wilcoxon"
cat("Random-reader random-case analysis")
cat("\nof Hupse Karssemeijer radiologist data:\n")
stats <- DiffFomAnal2007Hillis53 (
  dataset09,
  FOM = FOM 
)
cat(str(stats))
