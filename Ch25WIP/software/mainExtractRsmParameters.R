rm(list = ls())

stop("fix or delete me")

fileName <- "NICO" 
frocData <- loadDataFile(fileName)
retFileName <- paste0("ANALYZED/", "saveRetRoc", fileName)

load(retFileName)

S <- length(retSmRoc)
C <- length(retSmRoc)
for (j in 1:S){
  mu <- as.numeric(retSmRoc[[j]]$mu$mu)
  lambdaP <- as.numeric(retSmRoc[[j]]$lambdaP$lambdaP)
  nuP <- as.numeric(retSmRoc[[j]]$nuP$nuP)
  S[j] <- nuP * exp(-lambdaP)
  C[j] <- mu
  cat("rdr =", j, ", mu = ", mu, ", lambdaP= ", lambdaP, ", nuP = ", nuP, "\n")
}
tauSC <- cor(as.vector(S), as.vector(C), method = "kendall")
cat("tau S vs C = ", tauSC, "\n")
 
