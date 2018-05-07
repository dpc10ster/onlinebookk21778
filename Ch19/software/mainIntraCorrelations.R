# mainIntraCorrelations.R
# Fig. 19.8d,e,f,g,h,i
rm(list = ls())

library(foreach)
library(RJafroc)
library(doRNG)
library(doParallel)
library(ggplot2)

type <- "pearson";showPlot <- FALSE;showResultOneDataset <- FALSE
fileNames <-  c("TONY", "VD", "FR", 
                "FED", "JT", "MAG", 
                "OPT", "PEN", "NICO",
                "RUS", "DOB1", "DOB2", 
                "DOB3", "FZR")

avgS <- array(dim = length(fileNames))
avgC <-avgS;avgA <- avgS;RhoSC <- avgS;RhoSA <- avgS;RhoAC <- avgS
clusterParms <- list()
for (f in 1:length(fileNames)){
  if (showResultOneDataset) if (f != 11) next
  fileName <- fileNames[f]
  retFileName <- paste0("allResults", fileName) 
  sysAnalFileName <- system.file("ANALYZED/RSM6", retFileName, package = "RJafroc")
  if (file.exists(sysAnalFileName)){
    load(sysAnalFileName)
    I <- allResults[[1]]$I
    J <- allResults[[1]]$J
    S <- array(dim = c(I, J));C <- S;A <- S
    AllResIndx <- 0
    for (i in 1:I){
      for (j in 1:J){
        AllResIndx <- AllResIndx + 1
        mu <- allResults[[AllResIndx]]$retRsm$mu
        lambdaP <- allResults[[AllResIndx]]$retRsm$lambdaP
        nuP <- allResults[[AllResIndx]]$retRsm$nuP
        S[i, j] <- nuP * exp(-lambdaP)
        C[i, j] <- pnorm(mu/sqrt(2))
        A[i, j] <- allResults[[AllResIndx]]$retRsm$AUC
      }
    }
    
    avgS[f] <- mean(S)
    avgC[f] <- mean(C)
    avgA[f] <- mean(A)
    
    RhoSC[f] <- cor(as.vector(S), as.vector(C), method = type)
    RhoAC[f] <- cor(as.vector(A), as.vector(C), method = type)
    RhoSA[f] <- cor(as.vector(S), as.vector(A), method = type)
    cat("f = ", f,
        " S[f] =", avgS[f], 
        " C[f] =", avgC[f], 
        " A[f] =", avgA[f],
        " RhoSC[f] =", RhoSC[f], 
        " RhoSA[f] =", RhoSA[f], 
        " RhoAC[f] =", RhoAC[f],
        "\n")
    
    df <- data.frame(S = as.vector(S), C = as.vector(C))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = S, y = C)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x, size = 2) +
      geom_point(size = 5) +
      labs(title = ij) + 
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
      scale_color_manual(values = "black") + 
      theme(axis.title.y = element_text(size = 25,face="bold"),
            axis.title.x = element_text(size = 30,face="bold")) 
    if (showPlot) print(p)
    if (showResultOneDataset) {
      m <- lm(S ~ C, df);
      cat(fileNames[f], " S vs. C slope = ", as.numeric(m$coefficients[2]), ", r2 = ", summary(m)$r.squared,"\n")
    }
    
    df <- data.frame(S = as.vector(S), A = as.vector(A))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = S, y = A)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x, size = 2) +
      geom_point(size = 5) +
      labs(title = ij) + 
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
      scale_color_manual(values = "black") + 
      theme(axis.title.y = element_text(size = 25,face="bold"),
            axis.title.x = element_text(size = 30,face="bold")) 
    if (showPlot) print(p)
    if (showResultOneDataset) {
      m <- lm(S ~ A, df);
      cat(fileNames[f], " S vs. A slope = ", as.numeric(m$coefficients[2]), ", r2 = ", summary(m)$r.squared,"\n")
    }
    
    df <- data.frame(C = as.vector(C), A = as.vector(A))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = C, y = A)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x, size = 2) +
      geom_point(size = 5) +
      labs(title = ij) + 
      theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
      scale_color_manual(values = "black") + 
      theme(axis.title.y = element_text(size = 25,face="bold"),
            axis.title.x = element_text(size = 30,face="bold")) 
    if (showPlot) print(p)
    if (showResultOneDataset) {
      m <- lm(C ~ A, df);
      cat(fileNames[f], " C vs. A slope = ", as.numeric(m$coefficients[2]), ", r2 = ", summary(m)$r.squared,"\n")
    }
    
    clusterParms <- c(clusterParms, list(list(S = S, C = C, A = A)))
  }else{
    stop("Results file does not exist.")
  }
}
if (showResultOneDataset) stop("stop for one dataset analysis")

cat(
  "Avg S =", mean(avgS),"\n", 
  "Avg C =", mean(avgC),"\n", 
  "Avg A =", mean(avgA),"\n",
  "Avg rhoSC =", mean(RhoSC),"\n", 
  "Avg rhoSA =", mean(RhoSA),"\n", 
  "Avg rhoAC =", mean(RhoAC),"\n",
  "\n")
names(clusterParms) <- fileNames

if (!file.exists("IntraCorrelationBootstrapResults")){
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  B <- 200;seed <- 1
  rho <- foreach (b = 1:B, .options.RNG = seed, .combine = "rbind", .packages = "RJafroc") %dorng% {
    rhoSCb <- rep(NA, length(fileNames));rhoACb <- rhoSCb;rhoSAb <- rhoSCb
    Sb1 <- array(dim = c(length(fileNames)));Cb1 <- Sb1;Ab1 <- Sb1
    for (f in 1:length(fileNames)){
      fileName <- fileNames[f]
      retFileName <- paste0("allResults", fileName) 
      sysAnalFileName <- system.file("ANALYZED/RSM6", retFileName, package = "RJafroc")
      if (file.exists(sysAnalFileName)){
        load(sysAnalFileName)
        I <- length(clusterParms[[fileNames[f]]]$S[,1])
        J <- length(clusterParms[[fileNames[f]]]$S[1,])
        Sb <- array(dim = c(I,J,length(fileNames)));Cb <- Sb;Ab <- Sb
        
        jBs <- ceiling(runif(J) * J) # bootstrap readers
        
        Sb[,,f] <- clusterParms[[fileNames[f]]]$S[ , jBs]
        Cb[,,f] <- clusterParms[[fileNames[f]]]$C[ , jBs]
        Ab[,,f] <- clusterParms[[fileNames[f]]]$A[ , jBs]
        Sb1[f] <- mean(Sb[,,f]);Cb1[f] <- mean(Cb[,,f]);Ab1[f] <- mean(Ab[,,f])
        rhoSCb[f] <- cor(as.vector(Sb[,,f]), as.vector(Cb[,,f]), method = type)
        rhoACb[f] <- cor(as.vector(Ab[,,f]), as.vector(Cb[,,f]), method = type)
        rhoSAb[f] <- cor(as.vector(Sb[,,f]), as.vector(Ab[,,f]), method = type)
      }else{
        stop("Results file does not exist.")
      }
    }
    c(mean(rhoSCb), mean(rhoACb), mean(rhoSAb), mean(Sb1), mean(Cb1), mean(Ab1))
  }
  stopCluster(cl)
  save(bootStrapResults, file = "IntraCorrelationBootstrapResults")
} else load(file = "IntraCorrelationBootstrapResults")

rhoSCb <- data.frame(value = bootStrapResults[ , 1])
rhoACb <- data.frame(value = bootStrapResults[ , 2])
rhoSAb <- data.frame(value = bootStrapResults[ , 3])
Sb <- data.frame(value = bootStrapResults[ , 4])
Cb <- data.frame(value = bootStrapResults[ , 5])
Ab <- data.frame(value = bootStrapResults[ , 6])

histogram <- ggplot(Sb, aes(x = value)) + 
  geom_histogram(color = "white") + 
  xlab("S") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
  scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold")) 
print(histogram)
ciS <- quantile(Sb$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of S is", paste(ciS, collapse = ", "), "\n")

histogram <- ggplot(Cb, aes(x = value)) + 
  geom_histogram(color = "white") + 
  xlab("C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
  scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(histogram)
ciC <- quantile(Cb$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of C is", paste(ciC, collapse = ", "), "\n")

histogram <- ggplot(Ab, aes(x = value)) + 
  geom_histogram(color = "white") + 
  xlab("A") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
  scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(histogram)
ciA <- quantile(Ab$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of A is", paste(ciA, collapse = ", "), "\n")

histogram <- ggplot(rhoSCb, aes(x = value)) + 
  geom_histogram(color = "white") + 
  xlab("rho(S,C)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
  scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(histogram)
ciSC <- quantile(rhoSCb$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of the correlation between S and C is", paste(ciSC, collapse = ", "), "\n")

histogram <- ggplot(rhoACb, aes(x = value)) + 
  geom_histogram(color = "white") + 
  xlab("rho(A,C)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
  scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(histogram)
ciAC <- quantile(rhoACb$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of the correlation between A and C is", paste(ciAC, collapse = ", "), "\n")

histogram <- ggplot(rhoSAb, aes(x = value)) + 
  geom_histogram(color = "white") + 
  xlab("rho(S,A)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
  scale_color_manual(values = "black") + 
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"))
print(histogram)
ciSA <- quantile(rhoSAb$value, c(0.025, 0.975), type = 1)
cat("The empirical 95% CI of the correlation between S and A is", paste(ciSA, collapse = ", "), "\n")
