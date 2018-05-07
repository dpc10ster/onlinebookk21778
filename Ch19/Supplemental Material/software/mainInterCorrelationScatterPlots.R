rm(list = ls()) # mainInterCorrelationPlots.R #LOCK line numbers from now on
library(RJafroc)
library(ggplot2)
source("loadDataFile.R")

pathName <- "../../../06 E Online Appendices/E24 Datasets/"
type <- "pearson"
fileNames <-  c("TONY", "VD", "FR", "FED", "JT", "MAG", "OPT", "PEN", "NICO", 
                "RUS", "DOB1", "DOB2", "DOB3", "FZR")
for (f in 14:length(fileNames)){
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
    df <- data.frame(aucPro = as.vector(aucPro), aucRsm = as.vector(aucRsm))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = aucRsm, y = aucPro)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ 0 + x) +
      geom_point() + 
      labs(title = ij) 
    print(p)
    m <- lm(aucPro ~ 0 + aucRsm, df);
    cat(fileNames[f], " PRO vs. RSM AUCs slope = ", coef(m), ", r2 = ", summary(m)$r.squared,"\n")
    
    df <- data.frame(aucCbm = as.vector(aucCbm), aucRsm = as.vector(aucRsm))
    ij <- paste0("D", f, ", I = ", I, ", J = ", J)
    p <- ggplot(data = df, aes(x = aucRsm, y = aucCbm)) +
      geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ 0 + x) +
      geom_point() + 
      labs(title = ij) 
    print(p)
    m <- lm(aucCbm ~ 0 + aucRsm, df);
    cat(fileNames[f], " CBM vs. RSM AUCs slope = ", coef(m), ", r2 = ", summary(m)$r.squared,"\n")
    next
  }else{
    stop("Results file does not exist. Must analyze all datasets before running this.")
  }
}
