rm(list = ls()) #mainCompareFomOrderings.R
library(RJafroc)

stop("fix or delete me")

fileName <-  c("FED", "TONY", "VD", "JT", "FR", "MAG", "OPT", "RUS", "PEN", "VT1", "VT2", "FZR", "NICO")
fileName <- "FED" 
reportFileNamewAfroc <- paste0(fileName,"wAfroc.xlsx")
reportFileNamewHrRoc <- paste0(fileName,"HrAuc.xlsx")
frocData <- loadDataFile(fileName)
I <- length(frocData$NL[,1,1,1]);J <- length(frocData$NL[1,,1,1])
# wAFROC plots
wAfrocFoms <- FigureOfMerit(frocData)
# ROC plots
rocData <- FROC2HrROC(frocData)
rocFoms <- FigureOfMerit(rocData, fom = "Wilcoxon")

for (i in 1:I)
{
  for (j in 1:J)
  {
    for (ip in 1:I)
    {
      for (jp in 1:J)
      {
        if(j != jp) next
        prdct <- (wAfrocFoms[i,j] - wAfrocFoms[ip,jp])*(rocFoms[i,j] - rocFoms[ip,jp])
        next
      }
      
    }   
  }
    
}