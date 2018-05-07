rm(list = ls()) # delete all existing variables
library(RJafroc)

datafile <- ReadDataFile("1.xlsx")


stop("temp")
retDbmFroc  <- DBMHAnalysis(frocData) 

retDbmIroc  <- DBMHAnalysis(frocData, fom = "HrAuc") 
 
retDbmSongA1  <- DBMHAnalysis(frocData, fom = "SongA1") 

retDbmSongA2  <- DBMHAnalysis(frocData, fom = "SongA2") ## very slow

plotM <- c(1)
plotR <- c(1)
EmpiricalOpCharac(data = frocData, trts = plotM, rdrs = plotR, opChType = "FROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
EmpiricalOpCharac(data = frocData, trts = plotMAvg, rdrs = plotRAvg, OpChType = "ROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
EmpiricalOpCharac(data = frocData, trts = plotMAvg, rdrs = plotRAvg, OpChType = "AFROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
EmpiricalOpCharac(data = frocData, trts = plotMAvg, rdrs = plotRAvg, OpChType = "FROC")


