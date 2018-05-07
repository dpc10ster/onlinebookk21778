#install.packages(pkgs = "~/Desktop/RJafroc_0.0.1.tar.gz", repo = NULL, type = "source")
rm(list = ls())
library(RJafroc)

# ROC example
retDbmRoc  <- DBMHAnalysis(rocData, fom = "Wilcoxon") 

# FROC examples
# following three examples are for ROC data inferred from FROC data using different methods
retDbmHrAuc  <- DBMHAnalysis(frocData, fom = "HrAuc") # highest rating inferred ROC

## Not run: 
#retDbmSongA1  <- DBMHAnalysis(frocData, fom = "SongA1")
#retDbmSongA2  <- DBMHAnalysis(frocData, fom = "SongA2")
## End(Not run)

## default JAFROC analysis, wJAFROC FOM is assumed 
retDbmwJafroc  <- DBMHAnalysis(frocData)

## wJAFROC1 FOM (only to be used if there are no non-diseased cases)
retDbmwJafroc1  <- DBMHAnalysis(frocData, fom = "wJAFROC1")

retDbmJafroc  <- DBMHAnalysis(frocData, fom = "JAFROC")

## JAFROC1 FOM (only to be used if there are no non-diseased cases)
retDbmJafroc1  <- DBMHAnalysis(frocData, fom = "JAFROC1")

# ROI example
retDbmRoi  <- DBMHAnalysis(roiData, fom = "ROI")
