rm(list=ls(all=TRUE)) # MainCrossedModalities.R
library(RJafroc)

InputFilenameBase <- "CrossedModalitiesDataFile"
InputFilename <- paste0(InputFilenameBase,".xlsx")
fom <- "wAFROC"

cat("***Data file is ", InputFilenameBase, '***\n')
cat("fom = ", fom, '\n')

dataset <- DfReadCrossedModalities(InputFilename)
res1 <- StSignificanceTestingCrossedModalities(
  dataset,avgIndx = 1, fom = fom, option = "RRFC")
res2 <- StSignificanceTestingCrossedModalities(
  dataset,avgIndx = 2, fom = fom, option = "RRFC")
cat("**** Averaging over reconstruction index  *****\n")
print(res1)
cat("**** Averaging over mAs index  *****\n")
print(res2)