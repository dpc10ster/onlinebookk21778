rm(list = ls());#mainSsOrhCheck.R
library(RJafroc)
alpha <- 0.05
fileName <- "VanDyke.lrc"
#fileName <- "Franken1.lrc"
rocData <- ReadDataFile(fileName, format = "MRMC")
KStar <- length(rocData$NL[1,1,,1])
effectSizeArr <- c(0.05, 0.05)
cat("fileName = ", fileName, "\n")
for (e in 1:length(effectSizeArr))
{
  effectSize <- effectSizeArr[e]
  if (fileName == "Franken1.lrc") 
  {
    stop("temp")
    Cov1 <- 0.000351859
    Cov2 <- 0.000346505
    Cov3 <- 0.000221453
    Var <- 0.001393652
    msTR <- 0.000622731
    VarTR <- msTR - Var + Cov1 + max(Cov2-Cov3,0)
    VarTR <- max(VarTR,0)
  } else if (fileName == "VanDyke.lrc") {
    Cov1 <- 0.000351859
    Cov2 <- 0.000346505
    Cov3 <- 0.000221453
    Var <- 0.001393652
    msTR <- 0.000623
    VarTR <- msTR - Var + Cov1 + max(Cov2-Cov3,0)
    VarTR <- max(VarTR,0)
  }
  cat("VarTR = ", VarTR, "Var = ", Var, 
      ", Cov1 = ", Cov1, ", Cov2 = ", Cov2, ", Cov3  = ", Cov3, 
      ", effectSize = ", effectSize, "\n")
  
  if (fileName == "Franken1.lrc") 
  {
    JArr <- c(5)
    if (e == 1) KArr <- c(526, 294, 526) else KArr <- c(190, 107, 190)
  } else if (fileName == "VanDyke.lrc") {
    JArr <- c(5)
    if (e == 1) KArr <- c(2000, 526, 2000) else KArr <- c(266, 191, 266)
  }
  optionArr <- c("RRRC", "FRRC", "RRFC")
  for (j in 1:(length(JArr))) {
    for (k in 1:(length(KArr))) {
      J <- JArr[j]
      K <- KArr[k]
      J <- 8
      K <- 240
      DeltaNum <- J*effectSize^2/2
      cat("J = ", J, ", K = ", K, "\n")
      for (o in 1:(length(optionArr))) {
        Var1 <- Var;Cov11 <- Cov1;Cov21 <- Cov2;Cov31 <- Cov3;VarTR1 <- VarTR 
        option <- optionArr[o]
        if (option == "RRRC") {
          DeltaDenom <- max(VarTR,0)+(KStar/K)*(Var-Cov1+(J-1)*max(Cov2-Cov3,0))
          ncp <- DeltaNum/DeltaDenom    
          ddfNumSqBrkt <- (KStar/K)*(Var-Cov1+(J-1)*max(Cov2-Cov3,0))
          ddfDenomSqBrkt <- (KStar/K)*(Var-Cov1-max(Cov2-Cov3,0))
          ddf <- (J-1)*((VarTR+ddfNumSqBrkt)/(VarTR+ddfDenomSqBrkt))^2
        } else if (option == "FRRC") {
          DeltaDenom <- (KStar/K)*(Var-Cov1+(J-1)*max(Cov2-Cov3,0))
          ncp <- DeltaNum/DeltaDenom
          ddf <- (K-1)
        } else if (option == "RRFC") {
          DeltaDenom <- max(VarTR,0)+(KStar/K)*(Var-Cov1-max(Cov2-Cov3,0))
          ncp <- DeltaNum/DeltaDenom
          ddf <- (J-1)
        } else stop ("Bad value for option")
        FCrit <- qf(1 - alpha, 1, ddf)
        Power <- pf(FCrit, 1, ddf, ncp = ncp, FALSE)
        cat("option = ", option, ", FCrit = ", FCrit, "ddf = ", ddf, ", ncp = ", ncp, ", power = ", Power, "\n")
        next
      }
    }
  }
}