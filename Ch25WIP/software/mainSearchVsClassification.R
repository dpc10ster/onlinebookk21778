rm(list = ls()) # mainSearchVsClassification.R
library(foreach)
library(RJafroc)
library(doRNG)
library(doParallel)

stop("fix or delte me")
fileName <-  c("TONY", "VD", "FR", "FED", "JT", "MAG", "OPT", "PEN", "NICO", 
               "RUS", "DOB1", "DOB2", "DOB3", "FZR")
ret1 <- NULL
for (f in 1:length(fileName)){
  retFileName <- paste0("ANALYZED/URSM/", "saveRetRoc", fileName[f])
  if (file.exists(retFileName)){
    load(retFileName)
    frocData <- loadDataFile(fileName[f])
    I <- length(frocData$modalityID)
    J <- length(frocData$readerID)
    S <- array(dim = c(I, J))
    C <- S
    A <- S
    s <- 1
    for (i in 1:I){
      for (j in 1:J){
        mu <- as.numeric(retSmRoc[[s]]$mu$mu)
        lambdaP <- as.numeric(retSmRoc[[s]]$lambdaP$lambdaP)
        nuP <- as.numeric(retSmRoc[[s]]$nuP$nuP)
        
        S[i, j] <- nuP * exp(-lambdaP)
        C[i, j] <- mu
        A[i, j] <- as.numeric(retSmRoc[[s]]$AUC$AUC)
        s <- s + 1
      }
    }
    tauSC <- cor(as.vector(S), as.vector(C), method = "kendall")
    tauAC <- cor(as.vector(A), as.vector(C), method = "kendall")
    tauSA <- cor(as.vector(S), as.vector(A), method = "kendall")
    cat("Dataset:", fileName[f], ", I =", I, ", J =", J, 
        ", Avg S =", mean(S), ", Avg C =", mean(C), ", Avg A =", mean(A),
        ", tau(SvsC) =", tauSC, ", tau(AC) =", tauAC, ", tau(SA) =", tauSA, "\n")
    ret1 <- c(ret1, list(list(S = S, C = C, A = A)))
  }else{
    stop("Results file does not exist.")
  }
}
names(ret1) <- fileName

cl <- makeCluster(detectCores())
registerDoParallel(cl)
B <- 200
seed <- 1
tauSCBs <- rep(NA, B)
tauACBs <- tauSCBs
tauSABs <- tauACBs
tau <- foreach (b = 1:B, .options.RNG = seed, .combine = "rbind", .packages = "RJafroc") %dorng% {
  tauSCTemp <- rep(NA, length(fileName))
  tauACTemp <- tauSCTemp
  tauSATemp <- tauACTemp
  for (f in 1:length(fileName)){
    retFileName <- paste0("ANALYZED/URSM/saveRetRoc", fileName[f])
    if (file.exists(retFileName)){
      load(retFileName)
      frocData <- loadDataFile(fileName[f])
      I <- length(frocData$modalityID)
      J <- length(frocData$readerID)
      selectJ <- ceiling(runif(J) * J)
      Sb <- ret1[[fileName[f]]]$S[ , selectJ]
      Cb <- ret1[[fileName[f]]]$C[ , selectJ]
      Ab <- ret1[[fileName[f]]]$A[ , selectJ]
      tauSCTemp[f] <- cor(as.vector(Sb), as.vector(Cb), method = "kendall")
      tauACTemp[f] <- cor(as.vector(Ab), as.vector(Cb), method = "kendall")
      tauSATemp[f] <- cor(as.vector(Sb), as.vector(Ab), method = "kendall")
    }else{
      stop("Results file does not exist.")
    }
  }
  c(mean(tauSCTemp), mean(tauACTemp), mean(tauSATemp))
}
stopCluster(cl)
tauSCBs <- tau[ , 1]
tauACBs <- tau[ , 2]
tauSABs <- tau[ , 3]

hist(tauSCBs)
ciSC <- quantile(tauSCBs, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of the correlation between S and C is", ciSC, "\n")
