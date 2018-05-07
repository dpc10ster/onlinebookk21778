rm(list = ls()) #mainSsDbmhOpen.R
library(RJafroc)
alpha <- 0.05;cat("alpha = ", alpha, "\n")
fileName <- "VanDyke.lrc"
#fileName <- "Franken1.lrc"
rocData <- DfReadDataFile(
  fileName, 
  format = "MRMC")
retDbm <- StSignificanceTesting(
  dataset = rocData, 
  fom = "Wilcoxon", 
  method = "DBMH") 
varYTR <- retDbm$varComp$varComp[3]
varYTC <- retDbm$varComp$varComp[4]
varYEps <- retDbm$varComp$varComp[6]
effectSize <- retDbm$ciDiffTrtRRRC$Estimate
cat("effect size = ", effectSize, "\n")
#RRRC
J <- 10; K <- 163
ncp <- (0.5*J*K*(effectSize)^2)/
  (K*varYTR+max(J*varYTC,0)+varYEps)
MS <- UtilMeanSquares(rocData, 
                     fom = "Wilcoxon", 
                     method = "DBMH")
ddf <- (MS$msTR+max(MS$msTC-MS$msTRC,0))^2/
  (MS$msTR^2)*(J-1)
FCrit <- qf(1 - alpha, 1, ddf)
Power <- 1-pf(FCrit, 1, ddf, ncp = ncp)
cat("J =", J, "\nK =", K,"\nFCrit =", 
    FCrit, "\nddf =", ddf, "\nncp =", 
    ncp, "\nRRRC power =", Power, "\n")

#FRRC
J <- 10; K <- 133
ncp <- (0.5*J*K*(effectSize)^2)/
  (max(J*varYTC,0)+varYEps)
ddf <- (K-1)
FCrit <- qf(1 - alpha, 1, ddf)
Power <- 1-pf(FCrit, 1, ddf, ncp = ncp)
cat("J =", J, "\nK =", K,"\nFCrit =", FCrit, 
    "\nddf =", ddf, "\nncp =", 
    ncp, "\nFRRC power =", Power, "\n")

#RRFC
J <- 10; K <- 53
ncp <- (0.5*J*K*(effectSize)^2)/(K*varYTR+varYEps)
ddf <- (J-1)
FCrit <- qf(1 - alpha, 1, ddf)
Power <- 1-pf(FCrit, 1, ddf, ncp = ncp)
cat("J =", J, "\nK =", K,"\nFCrit =", FCrit, 
    "\nddf =", ddf, "\nncp =", 
    ncp, "\nRRFC power =", Power, "\n")
