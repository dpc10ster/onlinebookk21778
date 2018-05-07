rm(list = ls()) # mainBsRSMVsOthers.R # compares RSM to PROPROC and CBM # LOCK line numbers
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
clParms <- NULL # parameters passed to cluster
avgAucPro <- array(dim = length(fileNames))
avgAucRsm <- array(dim = length(fileNames))
avgAucCbm <- array(dim = length(fileNames))
avgSlopeProRsm <- array(dim = length(fileNames)) 
avgR2ProRsm <- array(dim = length(fileNames)) 
avgSlopeCbmRsm <- array(dim = length(fileNames)) 
avgR2CbmRsm <- array(dim = length(fileNames)) 
rhoMuRsmMuCbm <- array(dim = length(fileNames))
rhoNupRsmAlphaCbm <- array(dim = length(fileNames))
for (f in 1:length(fileNames)){
  retFileName <- paste0("ANALYZED/RSM-PROPROC-CBM/", "saveRetRoc", fileNames[f])
  if (file.exists(retFileName)){
    load(retFileName)
    frocData <- loadDataFile(fileNames[f], pathName)
    I <- length(frocData$modalityID)
    J <- length(frocData$readerID)
    aucRsm <- array(dim = c(I, J));aucCbm <- array(dim = c(I, J));aucPro <- array(dim = c(I, J))
    muRsm <- array(dim = c(I, J));muCbm <- array(dim = c(I, J))
    nupRsm <- array(dim = c(I, J));alphaCbm <- array(dim = c(I, J))
    s <- 1
    for (i in 1:I){
      for (j in 1:J){
        aucRsm[i, j] <- as.numeric(retSmRoc[[s]]$AUC$AUC)
        aucPro[i, j] <- as.numeric(retSmRoc[[s]]$aucProp)
        aucCbm[i, j] <- retSmRoc[[s]]$aucCbm
        muRsm[i, j] <- as.numeric(retSmRoc[[s]]$mu$mu)
        nupRsm[i, j] <- as.numeric(retSmRoc[[s]]$nuP$nuP)
        alphaCbm[i, j] <- retSmRoc[[s]]$alpha
        muCbm[i, j] <- retSmRoc[[s]]$muCbm
        s <- s + 1
      }
    }

    avgAucRsm[f] <- mean(aucRsm);avgAucPro[f] <- mean(aucPro);avgAucCbm[f] <- mean(aucCbm)
    rhoMuRsmMuCbm[f] <- cor(as.vector(muRsm), as.vector(muCbm),method = "pearson")
    rhoNupRsmAlphaCbm[f]<- cor(as.vector(nupRsm), as.vector(alphaCbm),method = "pearson")
    cat("rhoMuRsmMuCbm[f]=",rhoMuRsmMuCbm[f],
        ", rhoNupRsmAlphaCbm[f]=",rhoNupRsmAlphaCbm[f],"\n")
    
    df <- data.frame(aucRsm = as.vector(aucRsm), aucPro = as.vector(aucPro))
    m <- lm(aucPro ~ 0 + aucRsm, data = df)
    avgSlopeProRsm[f] <- coef(m)
    avgR2ProRsm[f] <- summary(m)$r.squared
    
    df <- data.frame(aucRsm = as.vector(aucRsm), aucCbm = as.vector(aucCbm))
    m <- lm(aucCbm ~ 0 + aucRsm, data = df)
    avgSlopeCbmRsm[f] <- coef(m)
    avgR2CbmRsm[f] <- summary(m)$r.squared
    
    clParms <- c(clParms, list(list(aucRsm = aucRsm, aucCbm = aucCbm, aucPro = aucPro, 
                                    nupRsm = nupRsm, alphaCbm = alphaCbm,
                                    muRsm = muRsm, muCbm = muCbm)))
  }else{
    stop("Results file does not exist. Must analyze all datasets before running this.")
  }
}
cat(
  "  avg aucRsm =", mean(avgAucRsm),
  ", avg aucPro =", mean(avgAucPro),
  ", avg aucCbm =", mean(avgAucCbm),
  ", avgSlopeCbmRsm =", mean(avgSlopeCbmRsm),", avg R2CbmRsm =", mean(avgR2CbmRsm),
  ", avgSlopeProRsm =", mean(avgSlopeProRsm),", avg R2ProRsm =", mean(avgR2ProRsm),
  ", avg rhoMuRsmMuCbm =", mean(rhoMuRsmMuCbm),", avg rhoNupRsmAlphaCbm =", mean(rhoNupRsmAlphaCbm),
  "\n"
)

# boostrap cluster code follows 
names(clParms) <- fileNames
cl <- makeCluster(detectCores())
registerDoParallel(cl)
B <- 200
seed <- 1
bootStrapResults <- foreach (b = 1:B, .options.RNG = seed, .combine = "rbind", .packages = "RJafroc") %dorng% {
  avgSlopeCbmRsm <- rep(NA, length(fileNames));avgR2CbmRsm <- rep(NA, length(fileNames))
  avgSlopeProRsm <- rep(NA, length(fileNames));avgR2ProRsm <- rep(NA, length(fileNames))
  rhoMuRsmMuCbm <- rep(NA, length(fileNames));rhoNupRsmAlphaCbm <- rep(NA, length(fileNames))
  for (f in 1:length(fileNames)){
    retFileName <- paste0("ANALYZED/RSM-PROPROC-CBM/", "saveRetRoc", fileNames[f])
    if (file.exists(retFileName)){
      load(retFileName)
      frocData <- loadDataFile(fileNames[f], pathName)
      I <- length(frocData$modalityID);J <- length(frocData$readerID)
      
      jBs <- ceiling(runif(J) * J) # here is were we bootstap readers
      
      # constrained fit thru origin; aucPro vs. aucRsm 
      df1 <- data.frame(aucRsm = as.vector(clParms[[fileNames[f]]]$aucRsm[ , jBs]), 
                       aucPro = as.vector(clParms[[fileNames[f]]]$aucPro[, jBs]))
      m <- lm(aucPro ~ 0 + aucRsm, data = df1)
      avgSlopeProRsm[f] <- coef(m)
      avgR2ProRsm[f] <- summary(m)$r.squared
      
      # constrained fit thru origin; aucCbm vs. aucRsm 
      df2 <- data.frame(aucRsm = as.vector(clParms[[fileNames[f]]]$aucRsm[ , jBs]), 
                       aucCbm = as.vector(clParms[[fileNames[f]]]$aucCbm[, jBs]))
      m <- lm(aucCbm ~ 0 + aucRsm, data = df2)
      avgSlopeCbmRsm[f] <- coef(m)
      avgR2CbmRsm[f] <- summary(m)$r.squared

      # correlation between muRsm and muCbm
      df1 <- data.frame(muRsm = as.vector(clParms[[fileNames[f]]]$muRsm[ , jBs]), 
                        muCbm = as.vector(clParms[[fileNames[f]]]$muCbm[, jBs]))
      rhoMuRsmMuCbm[f] <- cor(df1$muRsm, df1$muCbm,method = "pearson")
      
      # correlation between nupRsm and alphaCbm
      df2 <- data.frame(nupRsm = as.vector(clParms[[fileNames[f]]]$nupRsm[ , jBs]), 
                        alphaCbm = as.vector(clParms[[fileNames[f]]]$alphaCbm[, jBs]))
      rhoNupRsmAlphaCbm[f] <- cor(df2$nupRsm, df2$alphaCbm,method = "pearson")
      
    }else{
      stop("Results file does not exist.")
    }
  }
  c(mean(avgSlopeProRsm), mean(avgR2ProRsm), 
    mean(avgSlopeCbmRsm), mean(avgR2CbmRsm), mean(rhoMuRsmMuCbm), mean(rhoNupRsmAlphaCbm))
}
stopCluster(cl)

avgSlopeProRsm <- data.frame(value = bootStrapResults[ , 1])
avgR2ProRsm <- data.frame(value = bootStrapResults[ , 2])
avgSlopeCbmRsm <- data.frame(value = bootStrapResults[ , 3])
avgR2CbmRsm <- data.frame(value = bootStrapResults[ , 4])
rhoMuRsmMuCbm <- data.frame(value = bootStrapResults[ , 5])
rhoNupRsmAlphaCbm <- data.frame(value = bootStrapResults[ , 6])

histogram <- ggplot(avgSlopeProRsm, aes(x = value)) + geom_histogram(color = "white") +
  xlab("avgSlopeProRsm")
print(histogram)
ciAvgSlopeProRsm <- quantile(avgSlopeProRsm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of avgSlopeProRsm is", paste(ciAvgSlopeProRsm, collapse = ", "), "\n")

histogram <- ggplot(avgSlopeCbmRsm, aes(x = value)) + geom_histogram(color = "white") +
  xlab("avgSlopeCbmRsm")
print(histogram)
ciAvgSlopeCbmRsm <- quantile(avgSlopeCbmRsm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of avgSlopeCbmRsm is", paste(ciAvgSlopeCbmRsm, collapse = ", "), "\n")

histogram <- ggplot(rhoMuRsmMuCbm, aes(x = value)) + geom_histogram(color = "white") +
  xlab("rhoMuRsmMuCbm")
print(histogram)
ciRhoMuRsmMuCbm <- quantile(rhoMuRsmMuCbm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of rhoMuRsmMuCbm is", paste(ciRhoMuRsmMuCbm, collapse = ", "), "\n")

histogram <- ggplot(rhoNupRsmAlphaCbm, aes(x = value)) + geom_histogram(color = "white") +
  xlab("rhoNupRsmAlphaCbm")
print(histogram)
ciRhoNupRsmAlphaCbm <- quantile(rhoNupRsmAlphaCbm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of rhoNupRsmAlphaCbm is", paste(ciRhoNupRsmAlphaCbm, collapse = ", "), "\n")

df <- data.frame(nuP = nuP, alpha = alpha)
p <- ggplot(data = df, aes(x = nuP, y = alpha)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  geom_point()
print(p)
m <- lm(alpha ~ nuP, df);
cat("a = ", coef(m)[1],", b = ", coef(m)[2],", r2 = ", summary(m)$r.squared,"\n")
