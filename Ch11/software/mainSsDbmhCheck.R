rm(list = ls()) #mainSsDbmhCheck.R
library(RJafroc)
alpha <- 0.05
fileName <- "VanDyke.lrc"
cat("fileName = ", fileName, "\n")
effectSize <- 0.05
trY <- 0.0003686 # from Hillis 2004
tcY <- 0.016372
epsY <- 0.068255
powTab <- SsPowerTable(alpha = alpha, effectSize = effectSize, desiredPower = 0.8,  
                       method = "DBMH", option = "RRRC", trY, tcY, epsY)
print(powTab) 
# compare value for 8 readers, i.e., 249, to Hillis 2004 paper, power = 0.789 for 240 cases
# with 9 more cases it is reasonable that power will be 0.8
#     numReaders numCases power
#          8      249      0.8

fileName <- "Franken1.lrc"
cat("fileName = ", fileName, "\n")
effectSize <- 0.03
trY <- -0.000759 # from Hillis 2004; FOM is not stated; likely binormal according to earlier website doc
tcY <- -0.001393
epsY <- 0.083643
powTab <- SsPowerTable(alpha = alpha, effectSize = effectSize, desiredPower = 0.8,  
                       method = "DBMH", option = "RRRC", trY, tcY, epsY)
print(powTab)
# compare values to Hillis 2004 paper, Table 3, exact match
