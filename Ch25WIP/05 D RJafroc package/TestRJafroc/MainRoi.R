install.packages(pkgs = "~/Desktop/RJafroc_0.0.2.tar.gz", repo = NULL, type = "source")
library(RJafroc)
rm(list = ls())


##########################################################################################
### ROI paradigm

retORRoi  <- ORHAnalysis(roiData, fom = "ROI")
OutputReport(dataset = roiData, method = "DBMH", fom = "ROI", showWarnings = "FALSE")
OutputReport(dataset = roiData, method = "ORH", fom = "ROI", showWarnings = "FALSE")

OutputReport(dataset = roiData, method = "DBMH", fom = "ROI", showWarnings = "FALSE")
OutputReport(dataset = roiData, method = "ORH", fom = "ROI", showWarnings = "FALSE")


