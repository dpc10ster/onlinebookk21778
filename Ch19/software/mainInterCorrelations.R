# mainInterCorrelations.R 
# compares RSM to PROPROC and CBM
# uncomment line 29 to get Fig. 19.6, a,b
# Fig. 19.6c,d 
# Fig. 19.7a,b,c,d 
rm(list = ls()) 
library(foreach)
library(RJafroc)
library(doRNG)
library(doParallel)
library(ggplot2)

type <- "pearson";showPlot <- FALSE
fileNames <-  c("TONY", "VD", "FR", 
                "FED", "JT", "MAG", 
                "OPT", "PEN", "NICO",
                "RUS", "DOB1", "DOB2", 
                "DOB3", "FZR")
options(digits = 6)
clusterParms <- list() # parameters passed to cluster
avgAucPro <- array(dim = length(fileNames))
avgAucCbm <- avgAucPro;
avgAucRsm <- avgAucPro
avgSlopeProRsm <- avgAucPro 
avgSlopeCbmRsm <- avgAucPro 
avgR2ProRsm <- avgAucPro;avgR2CbmRsm <- avgAucPro 
rhoMuRsmMuCbm <- avgAucPro;rhoNupRsmAlphaCbm <- avgAucPro
for (f in 1:length(fileNames)){
  #if (f != 7) next
  fileName <- fileNames[f]
  retFileName <- paste0("allResults", fileName) 
  sysAnalFileName <- system.file("ANALYZED/RSM6", retFileName, package = "RJafroc")
  if (file.exists(sysAnalFileName)){
    load(sysAnalFileName)
    I <- allResults[[1]]$I
    J <- allResults[[1]]$J
    aucRsm <- array(dim = c(I, J));aucCbm <- array(dim = c(I, J));aucPro <- array(dim = c(I, J))
    muRsm <- array(dim = c(I, J));muCbm <- array(dim = c(I, J))
    nupRsm <- array(dim = c(I, J));alphaCbm <- array(dim = c(I, J))
    AllResIndx <- 0
    for (i in 1:I){
      for (j in 1:J){
        AllResIndx <- AllResIndx + 1
        aucRsm[i, j] <- allResults[[AllResIndx]]$retRsm$AUC
        aucPro[i, j] <- allResults[[AllResIndx]]$aucProp
        aucCbm[i, j] <- allResults[[AllResIndx]]$retCbm$AUC
        muRsm[i, j] <- allResults[[AllResIndx]]$retRsm$mu
        nupRsm[i, j] <- allResults[[AllResIndx]]$retRsm$nuP
        muCbm[i, j] <- allResults[[AllResIndx]]$retCbm$mu
        alphaCbm[i, j] <- allResults[[AllResIndx]]$retCbm$alpha
      }
    }
    
    # constrained fit thru origin; aucPro vs. aucRsm 
    df <- data.frame(aucRsm = as.vector(aucRsm), aucPro = as.vector(aucPro))
    mProRsm <- lm(aucPro ~ 0 + aucRsm, data = df)
    avgSlopeProRsm[f] <- coef(mProRsm)
    avgR2ProRsm[f] <- summary(mProRsm)$r.squared
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = aucRsm, y = aucPro)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ 0 + x) +
      geom_point(size = 5) + 
      labs(title = ij) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
      theme(
        axis.text=element_text(size=10),
        axis.title=element_text(size=28,face="bold"))
    if (showPlot) print(p)

    # constrained fit thru origin; aucCbm vs. aucRsm 
    df <- data.frame(aucRsm = as.vector(aucRsm), aucCbm = as.vector(aucCbm))
    mCbmRsm <- lm(aucCbm ~ 0 + aucRsm, data = df)
    avgSlopeCbmRsm[f] <- coef(mCbmRsm)
    avgR2CbmRsm[f] <- summary(mCbmRsm)$r.squared
    p <- ggplot(data = df, aes(x = aucRsm, y = aucCbm)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ 0 + x) +
      geom_point(size = 5) + 
      labs(title = ij) +
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
      theme(
        axis.text=element_text(size=10),
        axis.title=element_text(size=28,face="bold"))
    if (showPlot) print(p)
    
    df <- data.frame(muCbm = as.vector(muCbm), muRsm = as.vector(muRsm))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = muRsm, y = muCbm)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ 0 + x, size = 2) +
      geom_point(size = 5) +
      labs(title = ij) + 
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
      scale_color_manual(values = "black") + 
      theme(axis.title.y = element_text(size = 25,face="bold"),
            axis.title.x = element_text(size = 30,face="bold")) 
    if (showPlot) print(p)
    
    df <- data.frame(alphaCbm = as.vector(alphaCbm), nupRsm = as.vector(nupRsm))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = nupRsm, y = alphaCbm)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ 0 + x, size = 2) +
      geom_point(size = 5) +
      labs(title = ij) + 
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
      scale_color_manual(values = "black") + 
      theme(axis.title.y = element_text(size = 25,face="bold"),
            axis.title.x = element_text(size = 30,face="bold")) 
    if (showPlot) print(p)
    
    avgAucPro[f] <- mean(aucPro)
    avgAucCbm[f] <- mean(aucCbm)
    avgAucRsm[f] <- mean(aucRsm)
    
    rhoMuRsmMuCbm[f] <- cor(as.vector(muRsm), as.vector(muCbm),method = "pearson")
    rhoNupRsmAlphaCbm[f]<- cor(as.vector(nupRsm), as.vector(alphaCbm),method = "pearson")
    
    cat(avgAucPro[f],avgAucCbm[f],avgAucRsm[f],
        avgSlopeProRsm[f],avgR2ProRsm[f],
        avgSlopeCbmRsm[f],avgR2CbmRsm[f],
        rhoMuRsmMuCbm[f],rhoNupRsmAlphaCbm[f],"\n")
    
    
    clusterParms <- c(clusterParms, list(list(aucPro = aucPro, aucCbm = aucCbm, aucRsm = aucRsm,  
                                              nupRsm = nupRsm, alphaCbm = alphaCbm,
                                              muRsm = muRsm, muCbm = muCbm)))
  }else{
    stop("Results file does not exist. Must analyze all datasets with mainRSM.R before running this.")
  }
}
cat("\n")
cat(
  mean(avgAucPro),mean(avgAucCbm),mean(avgAucRsm),
  mean(avgSlopeProRsm),mean(avgR2ProRsm),
  mean(avgSlopeCbmRsm),mean(avgR2CbmRsm),
  mean(rhoMuRsmMuCbm),mean(rhoNupRsmAlphaCbm),"\n"
)
if (!file.exists("InterCorrelationBootstrapResults")){
  # boostrap cluster code follows 
  names(clusterParms) <- fileNames
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  B <- 200
  seed <- 1
  bootStrapResults <- foreach (b = 1:B, .options.RNG = seed, .combine = "rbind", .packages = "RJafroc") %dorng% {
    slopeCbmRsm <- rep(NA, length(fileNames));avgR2CbmRsm <- slopeCbmRsm
    slopeProRsm <- slopeCbmRsm;avgR2ProRsm <- slopeCbmRsm
    rhoMuRsmMuCbm <- rep(NA, length(fileNames));rhoNupRsmAlphaCbm <- rep(NA, length(fileNames))
    AucRsm1 <- array(dim = c(length(fileNames)));AucPro1 <- AucRsm1;AucCbm1 <- AucRsm1
    for (f in 1:length(fileNames)){
      fileName <- fileNames[f]
      retFileName <- paste0("allResults", fileName) 
      sysAnalFileName <- system.file("ANALYZED/RSM6", retFileName, package = "RJafroc")
      if (file.exists(sysAnalFileName)){
        load(sysAnalFileName)
        I <- allResults[[1]]$I
        J <- allResults[[1]]$J
        AucRsm <- array(dim = c(I,J,length(fileNames)));AucPro <- AucRsm;AucCbm <- AucRsm
        
        jBs <- ceiling(runif(J) * J) # here is where we bootstap readers
        
        AucRsm[,,f]  <-  clusterParms[[fileNames[f]]]$aucRsm[ , jBs]
        AucPro[,,f]  <-  clusterParms[[fileNames[f]]]$aucPro[ , jBs]
        AucCbm[,,f]  <-  clusterParms[[fileNames[f]]]$aucCbm[ , jBs]
        
        AucRsm1[f] <- mean(AucRsm[,,f]);AucPro1[f] <- mean(AucPro[,,f]);AucCbm1[f] <- mean(AucCbm[,,f])
        
        # constrained fit thru origin; aucPro vs. aucRsm 
        df1 <- data.frame(aucRsm = as.vector(clusterParms[[fileNames[f]]]$aucRsm[ , jBs]), 
                          aucPro = as.vector(clusterParms[[fileNames[f]]]$aucPro[, jBs]))
        mProRsm <- lm(aucPro ~ 0 + aucRsm, data = df1)
        slopeProRsm[f] <- coef(mProRsm);avgR2ProRsm[f] <- summary(mProRsm)$r.squared
        
        # constrained fit thru origin; aucCbm vs. aucRsm 
        df2 <- data.frame(aucRsm = as.vector(clusterParms[[fileNames[f]]]$aucRsm[ , jBs]), 
                          aucCbm = as.vector(clusterParms[[fileNames[f]]]$aucCbm[, jBs]))
        mCbmRsm <- lm(aucCbm ~ 0 + aucRsm, data = df2)
        slopeCbmRsm[f] <- coef(mCbmRsm);avgR2CbmRsm[f] <- summary(mCbmRsm)$r.squared
        
        # correlation between muRsm and muCbm
        df1 <- data.frame(muRsm = as.vector(clusterParms[[fileNames[f]]]$muRsm[ , jBs]), 
                          muCbm = as.vector(clusterParms[[fileNames[f]]]$muCbm[, jBs]))
        rhoMuRsmMuCbm[f] <- cor(df1$muRsm, df1$muCbm,method = "pearson")
        
        # correlation between nupRsm and alphaCbm
        df2 <- data.frame(nupRsm = as.vector(clusterParms[[fileNames[f]]]$nupRsm[ , jBs]), 
                          alphaCbm = as.vector(clusterParms[[fileNames[f]]]$alphaCbm[, jBs]))
        rhoNupRsmAlphaCbm[f] <- cor(df2$nupRsm, df2$alphaCbm,method = "pearson")
        
      }else{
        stop("Results file does not exist.")
      }
    }
    c(mean(AucPro1),mean(AucCbm1),mean(AucRsm1),
      mean(slopeProRsm),mean(avgR2ProRsm),
      mean(slopeCbmRsm),mean(avgR2CbmRsm),
      mean(rhoMuRsmMuCbm),mean(rhoNupRsmAlphaCbm))
  }
  stopCluster(cl)
  save(bootStrapResults, file = "InterCorrelationBootstrapResults")
} else load(file = "InterCorrelationBootstrapResults")

aucPro <- data.frame(value = bootStrapResults[ , 1])
aucCbm <- data.frame(value = bootStrapResults[ , 2])
aucRsm <- data.frame(value = bootStrapResults[ , 3])
slopeProRsm <- data.frame(value = bootStrapResults[ ,4])
#avgR2ProRsm <- data.frame(value = bootStrapResults[ ,5])
slopeCbmRsm <- data.frame(value = bootStrapResults[ ,6])
#avgR2CbmRsm <- data.frame(value = bootStrapResults[ ,7])
rhoMuRsmMuCbm <- data.frame(value = bootStrapResults[ ,8])
rhoNupRsmAlphaCbm <- data.frame(value = bootStrapResults[ ,9])

ciAvgAucPro <- quantile(aucPro$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of avgAucPro is", paste(ciAvgAucPro, collapse = ", "), "\n")

ciAvgAucCbm <- quantile(aucCbm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of avgAucCbm is", paste(ciAvgAucCbm, collapse = ", "), "\n")

ciAvgAucRsm <- quantile(aucRsm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of avgAucRsm is", paste(ciAvgAucRsm, collapse = ", "), "\n")

histogram <- ggplot(slopeProRsm, aes(x = value)) + geom_histogram(color = "white") +
  xlab("slopeProRsm") + scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) 
if (showPlot) print(histogram)
cislopeProRsm <- quantile(slopeProRsm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of slopeProRsm is", paste(cislopeProRsm, collapse = ", "), "\n")

histogram <- ggplot(slopeCbmRsm, aes(x = value)) + geom_histogram(color = "white") +
  xlab("slopeCbmRsm") + scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) 
if (showPlot) print(histogram)
cislopeCbmRsm <- quantile(slopeCbmRsm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of slopeCbmRsm is", paste(cislopeCbmRsm, collapse = ", "), "\n")

histogram <- ggplot(rhoMuRsmMuCbm, aes(x = value)) + 
  geom_histogram(color = "white") +
  xlab("rhoMuRsmMuCbm") + scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) 
if (showPlot) print(histogram)
ciRhoMuRsmMuCbm <- quantile(rhoMuRsmMuCbm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of rhoMuRsmMuCbm is", paste(ciRhoMuRsmMuCbm, collapse = ", "), "\n")

histogram <- ggplot(rhoNupRsmAlphaCbm, aes(x = value)) + 
  geom_histogram(color = "white") +
  xlab("rhoNupRsmAlphaCbm") + scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) 
if (showPlot) print(histogram)
ciRhoNupRsmAlphaCbm <- quantile(rhoNupRsmAlphaCbm$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of rhoNupRsmAlphaCbm is", paste(ciRhoNupRsmAlphaCbm, collapse = ", "), "\n")

