rm(list = ls()) # mainAnalyzeSingleModality.R
library(RJafroc);options(digits = 4)

# included datasets
fileNames <-  c("TONY", "VD", "FR", "FED", 
                "JT", "MAG", "OPT", "PEN", 
                "NICO","RUS", "DOB1", "DOB2", 
                "DOB3", "FZR")

f <- 4
fileName <- fileNames[f] # i.e., FED data
# the datasets already exist as R objects in RJafroc
theData <- get(sprintf("dataset%02d", f))

singleModData <- DfExtractDataset(
  theData, trts = 1)
singleModDataRes <- 
  StSignificanceTestingSingleFixedFactor(
    singleModData, fom = "wAFROC")
print(singleModDataRes)
# format(round(as.numeric(singleModDataRes$fomStats$stdErr), 4), nsmall = 4)
# format(round(as.numeric(singleModDataRes$fomStats$stdErr), 4), nsmall = 4)
# format(round(as.numeric(singleModDataRes$diffFomStats$stdErr), 4), nsmall = 4)

