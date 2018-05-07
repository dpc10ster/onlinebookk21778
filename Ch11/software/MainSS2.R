rm(list = ls()) #mainSS2.R
library(RJafroc)
fileName <- "VanDyke.lrc"
#fileName <- "Franken1.lrc"
rocData <- ReadDataFile(fileName, format = "MRMC")
retDbm <- DBMHAnalysis(dataset = rocData, fom = "Wilcoxon")                     
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
powTab <- PowerTable(dataset = rocData, effectSize = effectSize, alpha = 0.05)
print(powTab)
#DBMH variance components
tr <- retDbm$varComp$varComp[3]
tc <- retDbm$varComp$varComp[4]
trc <- retDbm$varComp$varComp[6]
J <- 10
K <- 163
delta <- (0.5*J*K*(effectSize)^2)/(K*tr+max(J*tc,0)+trc)
ddf <- (J-1)*(K*tr+max(J*tc,0)+trc)^2/(K*tr+trc)^2
alpha <- 0.05
FCrit <- qf(1 - alpha, 1, ddf)
Power <- 1-pf(FCrit, 1, ddf, ncp = delta)
cat("power = ", Power, "\n")

