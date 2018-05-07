rm(list = ls());#mainSsOrhDetails.R
library(RJafroc)

fileName <- "VanDyke"
fileName <- "Franken"
cat("fileName = ", fileName, "\n")
if (fileName == "VanDyke") rocData <- dataset02 else rocData <- dataset03
KStar <- length(rocData$NL[1,1,,1])
JStar <- length(rocData$NL[1,,1,1])
retORH <- StSignificanceTesting(rocData,FOM = "Wilcoxon", method = "ORH")
d <- retORH$ciDiffTrtRRRC$Estimate
Cov1 <- retORH$varComp$varCov[3]
Cov2 <- retORH$varComp$varCov[4]
Cov3 <- retORH$varComp$varCov[5]
Var <- retORH$varComp$varCov[6]
VarTR <- retORH$varComp$varCov[2]
VarTR <- max(VarTR,0)
cat("VarTR = ", VarTR, "Var = ", Var, 
    ", Cov1 = ", Cov1, ", Cov2 = ", Cov2, ", Cov3  = ", Cov3, ", d = ", d, "\n")

alpha <- 0.05
JArr <- c(JStar, 2*JStar)
KArr <- c(KStar, 2*KStar)
optionArr <- c("RRRC", "FRRC", "RRFC")
for (j in 1:(length(JArr))) {
  for (k in 1:(length(KArr))) {
    J <- JArr[j]
    K <- KArr[k]
    DeltaNum <- J*d^2/2
    cat("J = ", J, ", K = ", K, "\n")
    for (o in 1:(length(optionArr))) {
      option <- optionArr[o]
      if (option == "RRRC") {
        DeltaDenom <- max(VarTR,0)+(KStar/K)*(Var-Cov1+(J-1)*max(Cov2-Cov3,0)) 
        Delta <- DeltaNum/DeltaDenom    
        ddfNumSqBrkt <- (KStar/K)*(Var-Cov1+(J-1)*max(Cov2-Cov3,0))
        ddfDenomSqBrkt <- (KStar/K)*(Var-Cov1-max(Cov2-Cov3,0))
        ddf <- (J-1)*((VarTR+ddfNumSqBrkt)/(VarTR+ddfDenomSqBrkt))^2
      } else if (option == "FRRC") {
        DeltaDenom <- (KStar/K)*(Var-Cov1+(J-1)*max(Cov2-Cov3,0))
        Delta <- DeltaNum/DeltaDenom    
        ddf <- (K-1)
      } else if (option == "RRFC") {
        DeltaDenom <- max(VarTR,0)+(KStar/K)*(Var-Cov1-max(Cov2-Cov3,0))
        Delta <- DeltaNum/DeltaDenom
        ddf <- (J-1)
      } else stop ("Bad value for option")
      power <- SsPowerGivenJK(
        J,K,
        effectSize = d,
        option = option,
        method = "ORH",
        cov1 = Cov1,
        cov2 = Cov2,
        cov3 = Cov3,
        varTR = VarTR,
        varEps = Var,
        KStar = KStar)
      if (o == 1) power <- power$powerRRRC
      if (o == 2) power <- power$powerFRRC
      if (o == 3) power <- power$powerRRFC
      FCrit <- qf(1 - alpha, 1, ddf)
      Power <- pf(FCrit, 1, ddf, ncp = Delta, FALSE)
      cat(option, "FCrit = ", FCrit, ", ddf = ", ddf, ", Delta = ", Delta, ", power = ", Power,
          ", Power RJafroc = ", power, "\n")
      next
    }
  }
}
