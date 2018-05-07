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
reportFileNamewAfroc <- paste0(fileName,"wAfroc.xlsx")
reportFileNamewHrRoc <- paste0(fileName,"HrAuc.xlsx")

# wAFROC plots
# plots for all treatment reader combinations
plotT <- c(4,5);plotR <- c(1,2,3,4) 
p <- PlotEmpiricalOperatingCharacteristics(theData,trts = plotT, rdrs = plotR, opChType = "wAFROC")$Plot
p <- p + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(0.7,0.15), legend.direction = "horizontal",
        legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(0.75, "inch")) 
p$layers[[1]]$aes_params$size <- 2 # line
#p$layers[[2]]$aes_params$size <- 5 # points
print(p)

plotT <- list(4,5) # create reader averaged curves in each modality
plotR <- list(c(1:4),c(1:4)) # create reader averaged curves in each modality
p <- PlotEmpiricalOperatingCharacteristics(theData,trts = plotT, rdrs = plotR, opChType = "wAFROC")$Plot
p <- p + 
  scale_color_manual(values = c("black","darkgrey")) +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(0.7,0.15), legend.direction = "horizontal",
        legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(0.75, "inch")) 
p$layers[[1]]$aes_params$size <- 2 # line
#p$layers[[2]]$aes_params$size <- 5 # points
print(p)

# ROC plots
# plots for all treatment reader combinations 
plotT <- c(4,5);plotR <- c(1,2,3,4)
rocData <- DfFroc2Roc(theData)
p <- PlotEmpiricalOperatingCharacteristics(rocData,trts = plotT, rdrs = plotR,opChType = "ROC")$Plot
p <- p + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(0.7,0.15), legend.direction = "horizontal",
        legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(0.75, "inch")) 
p$layers[[1]]$aes_params$size <- 2 # line
#p$layers[[2]]$aes_params$size <- 5 # points
print(p)

plotT <- list(4,5) # trick to create reader averaged curves
plotR <- list(c(1:4),c(1:4)) # trick to create reader averaged curve
p <- PlotEmpiricalOperatingCharacteristics(rocData,trts = plotT, rdrs = plotR,opChType = "ROC")$Plot
p <- p + 
  scale_color_manual(values = c("black","darkgrey")) +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(0.7,0.15), legend.direction = "horizontal",
        legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(0.75, "inch")) 
p$layers[[1]]$aes_params$size <- 2 # line
print(p)
