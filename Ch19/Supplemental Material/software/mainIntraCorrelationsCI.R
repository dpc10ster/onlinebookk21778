rm(list = ls()) # mainIntraCorrelationsCI.R # LOCK line numbers
library(foreach)
library(RJafroc)
library(doRNG)
library(doParallel)
library(ggplot2)
source("loadDataFile.R")

pathName <- "../../../06 E Online Appendices/E24 Datasets/"
type <- "pearson"
fileNames <-  c("TONY", "VD", "FR", "FED", "JT", "MAG", "OPT", "PEN", "NICO", 
               "RUS", "DOB1", "DOB2", "DOB3", "FZR")
avgS <- array(dim = length(fileNames))
avgAucC <-avgS;avgAucRsm <- avgS;RhoSC <- avgS;RhoSA <- avgS;RhoAC <- avgS
clParms <- NULL
for (f in 1:length(fileNames)){
  retFileName <- paste0("ANALYZED/RSM-PROPROC-CBM/", "saveRetRoc", fileNames[f])
  if (file.exists(retFileName)){
    load(retFileName)
    frocData <- loadDataFile(fileNames[f], pathName)
    I <- length(frocData$modalityID)
    J <- length(frocData$readerID)
    S <- array(dim = c(I, J));C <- S;A <- S
    s <- 1
    for (i in 1:I){
      for (j in 1:J){
        mu <- as.numeric(retSmRoc[[s]]$mu$mu)
        lambdaP <- as.numeric(retSmRoc[[s]]$lambdaP$lambdaP)
        nuP <- as.numeric(retSmRoc[[s]]$nuP$nuP)
        
        S[i, j] <- nuP * exp(-lambdaP)
        C[i, j] <- pnorm(mu/sqrt(2))
        A[i, j] <- as.numeric(retSmRoc[[s]]$AUC$AUC)
        s <- s + 1
      }
    }
    avgS[f] <- mean(S)
    avgAucC[f] <- mean(C)
    avgAucRsm[f] <- mean(A)
    rhoSC <- cor(as.vector(S), as.vector(C), method = type)
    rhoAC <- cor(as.vector(A), as.vector(C), method = type)
    rhoSA <- cor(as.vector(S), as.vector(A), method = type)
    RhoSC[f] <- rhoSC
    RhoSA[f] <- rhoSA
    RhoAC[f] <- rhoAC
    df <- data.frame(aucRsm = as.vector(A))
    clParms <- c(clParms, list(list(S = S, C = C, A = A)))
  }else{
    stop("Results file does not exist.")
  }
}
cat(
  "  Avg S =", mean(avgS), 
  ", Avg C =", mean(avgAucC), 
  ", Avg A =", mean(avgAucRsm),
  "  Avg rhoSC =", mean(RhoSC), 
  ", Avg rhoSA =", mean(RhoSA), 
  ", Avg rhoAC =", mean(RhoAC),
  "\n")
names(clParms) <- fileNames

cl <- makeCluster(detectCores())
registerDoParallel(cl)
B <- 200;seed <- 1
rho <- foreach (b = 1:B, .options.RNG = seed, .combine = "rbind", .packages = "RJafroc") %dorng% {
  rhoSCTemp <- rep(NA, length(fileNames));rhoACTemp <- rhoSCTemp;rhoSATemp <- rhoSCTemp
  Sb1 <- array(dim = c(length(fileNames)));Cb1 <- Sb1;Ab1 <- Sb1
  for (f in 1:length(fileNames)){
    retFileName <- paste0("ANALYZED/RSM-PROPROC-CBM/", "saveRetRoc", fileNames[f])
    if (file.exists(retFileName)){
      I <- length(clParms[[fileNames[f]]]$S[,1])
      J <- length(clParms[[fileNames[f]]]$S[1,])
      Sb <- array(dim = c(I,J,length(fileNames)));Cb <- Sb;Ab <- Sb
      
      jBs <- ceiling(runif(J) * J) # bootstrap readers
      
      Sb[,,f] <- clParms[[fileNames[f]]]$S[ , jBs]
      Cb[,,f] <- clParms[[fileNames[f]]]$C[ , jBs]
      Ab[,,f] <- clParms[[fileNames[f]]]$A[ , jBs]
      Sb1[f] <- mean(Sb[,,f]);Cb1[f] <- mean(Cb[,,f]);Ab1[f] <- mean(Ab[,,f])
      rhoSCTemp[f] <- cor(as.vector(Sb[,,f]), as.vector(Cb[,,f]), method = type)
      rhoACTemp[f] <- cor(as.vector(Ab[,,f]), as.vector(Cb[,,f]), method = type)
      rhoSATemp[f] <- cor(as.vector(Sb[,,f]), as.vector(Ab[,,f]), method = type)
    }else{
      stop("Results file does not exist.")
    }
  }
  c(mean(rhoSCTemp), mean(rhoACTemp), mean(rhoSATemp), mean(Sb1), mean(Cb1), mean(Ab1))
}
stopCluster(cl)

rhoSCBs <- data.frame(value = rho[ , 1])
rhoACBs <- data.frame(value = rho[ , 2])
rhoSABs <- data.frame(value = rho[ , 3])
S <- data.frame(value = rho[ , 4])
C <- data.frame(value = rho[ , 5])
A <- data.frame(value = rho[ , 6])

histogram <- ggplot(S, aes(x = value)) + geom_histogram(color = "white") + xlab("S")
print(histogram)
ciS <- quantile(S$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of S is", paste(ciS, collapse = ", "), "\n")

histogram <- ggplot(C, aes(x = value)) + geom_histogram(color = "white") + xlab("C")
print(histogram)
ciC <- quantile(C$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of C is", paste(ciC, collapse = ", "), "\n")

histogram <- ggplot(A, aes(x = value)) + geom_histogram(color = "white") + xlab("A")
print(histogram)
ciA <- quantile(A$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of A is", paste(ciA, collapse = ", "), "\n")

histogram <- ggplot(rhoSCBs, aes(x = value)) + geom_histogram(color = "white") + xlab("rho(S,C)")
print(histogram)
ciSC <- quantile(rhoSCBs$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of the correlation between S and C is", paste(ciSC, collapse = ", "), "\n")

histogram <- ggplot(rhoACBs, aes(x = value)) + geom_histogram(color = "white") + xlab("rho(A,C)")
print(histogram)
ciAC <- quantile(rhoACBs$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of the correlation between A and C is", paste(ciAC, collapse = ", "), "\n")

histogram <- ggplot(rhoSABs, aes(x = value)) + geom_histogram(color = "white") + xlab("rho(S,A)")
print(histogram)
ciSA <- quantile(rhoSABs$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of the correlation between S and A is", paste(ciSA, collapse = ", "), "\n")
