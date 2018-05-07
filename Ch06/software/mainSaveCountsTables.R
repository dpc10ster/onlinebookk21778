rm( list = ls())#mainSaveCountsTables.R
library(RJafroc)

K1 <- c(30,19,8,2,1) # this is the observed data!
K2 <- c(5,6,5,12,22) # this is the observed data!
# convert to linear vectors
NL <- NULL;LL <- NULL
for (i in 1:length(K1)) {
  NL <- c(NL,rep(i, K1[i]))
}
for (i in 1:length(K2)) {
  LL <- c(LL,rep(i, K2[i]))
}

rocData <- Df2RJafrocDataset(NL, LL)
plotRoc <- PlotEmpiricalOperatingCharacteristics(rocData, 1, 1)
print(plotRoc$Plot)
DfSaveDataFile(rocData,fileName = "Table 4.1.xlsx", format = "JAFROC" )
K1 <- c(30,19,8,7,5) # this is the cheated data!
K2 <- c(5,6,5,12,22) # this is the observed data!
NL <- NULL;LL <- NULL
for (i in 1:length(K1)) {
  NL <- c(NL,rep(i, K1[i]))
}
for (i in 1:length(K2)) {
  LL <- c(LL,rep(i, K2[i]))
}
rocData <- Df2RJafrocDataset(NL, LL)
plotRoc <- PlotEmpiricalOperatingCharacteristics(rocData, 1, 1)
print(plotRoc$Plot)
DfSaveDataFile(rocData,fileName = "Table 4.1Cheat.xlsx", format = "JAFROC" )
