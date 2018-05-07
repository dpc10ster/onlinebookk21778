rm( list = ls()) # mainWilcoxon.R
library(caTools)
source("Wilcoxon.R");source("RocOperatingPoints.R")
options(digits = 9)
RocCountsTable = array(dim = c(2,5))
RocCountsTable[1,]  <- c(30,19,8,2,1)
RocCountsTable[2,]  <- c(5,6,5,12,22)

#convert frequency table to array
zk1  <- rep(1:length(RocCountsTable[1,]),
            RocCountsTable[1,])
zk2  <- rep(1:length(RocCountsTable[2,]),
            RocCountsTable[2,])

w  <- Wilcoxon (zk1, zk2)
cat("The wilcoxon statistic is = ", w, "\n")
ret <- RocOperatingPoints(RocCountsTable[1,], 
                          RocCountsTable[2,])
FPF <- ret$FPF;FPF <- c(0,FPF,1)
TPF <- ret$TPF;TPF <- c(0,TPF,1)
AUC <- trapz(FPF,TPF) # trapezoidal integration
cat("direct integration yields AUC = ", AUC, "\n")



