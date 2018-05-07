rm(list = ls()) #MainBetweenReaderSd.R
library(RJafroc)
# dataset03 contains the Franken dataset
# dataset02 contains the VanDyke dataset
# dataRoc <- ReadDataFile(fileName, format = "MRMC") # this is not needed as dataset already exists in RJafroc
Foms <- UtilFigureOfMerit(dataset02, fom = "Wilcoxon")
var1 <- var(Foms[1,]);var2 <- var(Foms[2,])
cat("between-reader variance in modality 1 =", var1, 
    "\nbetween-reader variance in modality 2 =", var2, 
    "\navg. between-reader variance in both modalities =", (var1+var2)/2, "\n")