rm(list = ls()) #mainwAfrocPlots.R
library(RJafroc);library(ggplot2)

# included datasets
fileNames <-  c("TONY", "VD", "FR", "FED", 
                "JT", "MAG", "OPT", "PEN", 
                "NICO","RUS", "DOB1", "DOB2", 
                "DOB3", "FZR")

f <- 4
fileName <- fileNames[f] # i.e., FED data
# the datasets already exist as R objects in RJafroc
theData <- get(sprintf("dataset%02d", f))

# wAFROC plots
# plots for all treatment reader combinations
plotT <- c(4,5);plotR <- c(1,2,3,4) 
p <- PlotEmpiricalOperatingCharacteristics(
  theData,trts = plotT, rdrs = plotR, 
  opChType = "wAFROC")
print(p$Plot)

# create reader averaged curves in each modality
plotT <- list(4,5);plotR <- list(c(1:4),c(1:4))
p <- PlotEmpiricalOperatingCharacteristics(
  theData,trts = plotT, rdrs = plotR, 
  opChType = "wAFROC")
print(p$Plot)

# ROC plots
# plots for all treatment reader combinations 
plotT <- c(4,5);plotR <- c(1,2,3,4)
rocData <- DfFroc2Roc(theData)
p <- PlotEmpiricalOperatingCharacteristics(
  rocData,trts = plotT, rdrs = plotR,
  opChType = "ROC")
print(p$Plot)

plotT <- list(4,5);plotR <- list(c(1:4),c(1:4))
p <- PlotEmpiricalOperatingCharacteristics(
  rocData,trts = plotT, rdrs = plotR,
  opChType = "ROC")
print(p$Plot)
