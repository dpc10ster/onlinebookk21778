rm(list = ls()) #MainModalityEffect.R
library(RJafroc)
# dataset03 contains the Franken dataset
# dataset02 contains the VanDyke dataset
# dataRoc <- ReadDataFile(fileName, format = "MRMC") # this is not needed as dataset already exists in RJafroc
Foms <- UtilFigureOfMerit(dataset02, FOM = "Wilcoxon")
mean1 <- mean(Foms[1,]);mean2 <- mean(Foms[2,])
cat("reader-average FOM in modality 1 =", mean1, "\nreader-average FOM in modality 2 =", mean2, 
    "\neffect size, i.e., fom modality 1 minus modality 2 =", mean1-mean2, "\n")

