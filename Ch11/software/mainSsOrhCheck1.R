rm(list = ls());#mainSsOrhCheck1.R
library(RJafroc)
alpha <- 0.05
fileName <- "VanDyke.lrc"
rocData <- ReadDataFile(fileName, format = "MRMC")
KStar <- length(rocData$NL[1,1,,1])
effectSize <- 0.05
cat("fileName = ", fileName, "\n")
# following are proproc generated values from Hillis 2011 for Van Dyke Data
Cov1 <- 0.000351859
Cov2 <- 0.000346505
Cov3 <- 0.000221453
Var <- 0.001393652
msTR <- 0.000623
VarTR <- msTR - Var + Cov1 + max(Cov2-Cov3,0)
VarTR <- max(VarTR,0)
cat("VarTR = ", VarTR, "Var = ", Var, 
    ", Cov1 = ", Cov1, ", Cov2 = ", Cov2, ", Cov3  = ", Cov3, 
    ", effectSize = ", effectSize, "\n")

optionArr <- c("RRRC")#, "FRRC", "RRFC")
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
