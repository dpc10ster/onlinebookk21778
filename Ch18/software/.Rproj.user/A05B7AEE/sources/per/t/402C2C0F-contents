rm(list = ls()) #mainAnalyzewAFROC.R
library(RJafroc)

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
if (!file.exists(reportFileNamewAfroc)){
  UtilOutputReport(
    dataset = theData, 
    reportFormat = "xlsx", 
    reportFile = reportFileNamewAfroc, 
    renumber = "TRUE",
    method = "DBMH", fom = "wAFROC", 
    showWarnings = "FALSE")
}
if (!file.exists(reportFileNamewHrRoc)){
  UtilOutputReport(
    dataset = theData, 
    reportFormat = "xlsx", 
    reportFile = reportFileNamewHrRoc, 
    renumber = "TRUE",
    method = "DBMH", fom = "HrAuc", 
    showWarnings = "FALSE")
}