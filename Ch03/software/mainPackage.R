rm(list = ls()) #mainPackage.R
library(RJafroc)
fileName <- "VanDyke.lrc"
rocData <- ReadDataFile(fileName, format = "MRMC")
