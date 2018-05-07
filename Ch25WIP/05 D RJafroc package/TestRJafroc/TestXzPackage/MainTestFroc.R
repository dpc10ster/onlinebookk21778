install.packages(pkgs = "~/Desktop/RJafroc_0.1.1.tar.gz", repo = NULL, type = "source")

rm(list = ls())
library(RJafroc)

retVortexOrg <- ReadDataFile(fileName="../Datasets/JAFROC_input_3mm_to_20mm.xlsx")
fom <- FigureOfMerit(dataset = retVortexOrg)
diffFomOrg <- mean(fom[,c(3,5)])-mean(fom[,c(1,2,4)]) 
K <- length(retVortexOrg$NL[1,1,,1])
K2 <- length(retVortexOrg$LL[1,1,,1])

retVortex <- retVortexOrg
B <- 200
diffFom <- array(B)
for (b in 1:B){
  kb <- ceiling(runif(K) * K)
  kb2 <- ceiling(runif(K2) * K2)
  NLb <- retVortexOrg$NL[,,kb,]
  LLb <- retVortexOrg$LL[,,kb2,]
  retVortex$NL <- NLb
  retVortex$LL <- LLb
  fom <- FigureOfMerit(dataset = retVortex)
  diffFom[b] <- mean(fom[,c(3,5)])-mean(fom[,c(1,2,4)]) 
}

diffFomCI <- quantile(diffFom, c(0.025,0.975))

retDbmFroc  <- DBMHAnalysis(frocData) 
print(retDbmFroc)

retDbmIroc  <- DBMHAnalysis(frocData, fom = "HrAuc") 
print(retDbmIroc)
 
retDbmSongA1  <- DBMHAnalysis(frocData, fom = "SongA1") 
print(retDbmSongA1)

retDbmSongA2  <- DBMHAnalysis(frocData, fom = "SongA2") 
print(retDbmSongA2)

plotM <- c(1:2)
plotR <- c(1:4)
EmpiricalOpCharac(data = frocData, modalities = plotM, readers = plotR, OpChType = "ROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
EmpiricalOpCharac(data = frocData, modalities = plotMAvg, readers = plotRAvg, OpChType = "ROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
EmpiricalOpCharac(data = frocData, modalities = plotMAvg, readers = plotRAvg, OpChType = "AFROC")

plotMAvg <- list(1, 2)
plotRAvg <- list(c(1:4),c(1:4))
EmpiricalOpCharac(data = frocData, modalities = plotMAvg, readers = plotRAvg, OpChType = "FROC")


