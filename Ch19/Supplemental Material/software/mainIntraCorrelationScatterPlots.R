rm(list = ls()) # mainIntraCorrelationScatterPlots.R # LOCK line numbers
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
avgC <- array(dim = length(fileNames))
avgA <- array(dim = length(fileNames))
RhoSC <- array(dim = length(fileNames))
RhoSA <- array(dim = length(fileNames))
RhoAC <- array(dim = length(fileNames))
for (f in 14:length(fileNames)){
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
    
    df <- data.frame(S = as.vector(S), C = as.vector(C))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = S, y = C)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~  x) +
      geom_point() +
      labs(title = ij)
    print(p)
    # m <- lm(S ~ C, df);
    # cat(fileNames[f], " S vs. C slope = ", as.numeric(m$coefficients[2]), ", r2 = ", summary(m)$r.squared,"\n")

    df <- data.frame(S = as.vector(S), A = as.vector(A))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = S, y = A)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~  x) +
      geom_point() +
      labs(title = ij)
    print(p)
    #m <- lm(S ~ A, df);
    #cat(fileNames[f], " S vs. A slope = ", as.numeric(m$coefficients[2]), ", r2 = ", summary(m)$r.squared,"\n")

    df <- data.frame(C = as.vector(C), A = as.vector(A))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = C, y = A)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~  x) +
      geom_point() +
      labs(title = ij)
    print(p)
    # m <- lm(C ~ A, df);
    # cat(fileNames[f], " C vs. A slope = ", as.numeric(m$coefficients[2]), ", r2 = ", summary(m)$r.squared,"\n")
    
    avgS[f] <- mean(S)
    avgC[f] <- mean(C)
    avgA[f] <- mean(A)
    
    RhoSC[f] <- cor(as.vector(S), as.vector(C), method = type)
    RhoAC[f] <- cor(as.vector(A), as.vector(C), method = type)
    RhoSA[f] <- cor(as.vector(S), as.vector(A), method = type)
    cat(
      "Avg S[f] =", avgS[f], 
      ", Avg C[f] =", avgC[f], 
      ", Avg A[f] =", avgA[f],
      ", RhoSC[f] =", RhoSC[f], 
      ", RhoSA[f] =", RhoSA[f], 
      ", RhoAC[f] =", RhoAC[f],
      "\n")
  stop("comment this to see other datasets")  
  }else{
    stop("Results file does not exist; need to run mainRSM.R after running PROPROC on dataset.")
  }
}
