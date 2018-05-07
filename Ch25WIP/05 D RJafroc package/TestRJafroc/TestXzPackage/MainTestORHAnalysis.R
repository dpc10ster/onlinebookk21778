#install.packages(pkgs = "~/Desktop/RJafroc_0.0.1.tar.gz", repo = NULL, type = "source")
rm(list = ls())
library(RJafroc)

# ROC example
retORRoc  <- ORHAnalysis(rocData, fom = "Wilcoxon") 

# FROC examples
# following three examples are for ROC data inferred from FROC data using different methods
retORHrAuc  <- ORHAnalysis(frocData, fom = "HrAuc") # highest rating inferred ROC

## Not run: 
#retORSongA1  <- ORHAnalysis(frocData, fom = "SongA1")
#retORSongA2  <- ORHAnalysis(frocData, fom = "SongA2")
## End(Not run)

## default JAFROC analysis, wJAFROC FOM is assumed 
retORwJafroc  <- ORHAnalysis(frocData)

## wJAFROC1 FOM (only to be used if there are no non-diseased cases)
retORwJafroc1  <- ORHAnalysis(frocData, fom = "wJAFROC1")

retORJafroc  <- ORHAnalysis(frocData, fom = "JAFROC")

## JAFROC1 FOM (only to be used if there are no non-diseased cases)
retORJafroc1  <- ORHAnalysis(frocData, fom = "JAFROC1")

# ROI example
retORRoi  <- ORHAnalysis(roiData, fom = "ROI")
