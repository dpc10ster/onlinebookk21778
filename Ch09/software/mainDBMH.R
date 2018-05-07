rm(list = ls()) # mainDBMH.R # freeze lines
library(RJafroc)

alpha <- 0.05
#fileName <- "Franken1.lrc"
fileName <- "VanDyke.lrc"
rocData <- DfReadDataFile(fileName, format = "MRMC")
pseudoValues <- UtilPseudoValues(rocData, FOM = "Wilcoxon")
I <- dim(pseudoValues)[1]
J <- dim(pseudoValues)[2]
K <- dim(pseudoValues)[3]
retMS <- UtilMeanSquares(rocData, FOM = "Wilcoxon")
FOM <- UtilFigureOfMerit(rocData, FOM = "Wilcoxon")

cat("\nRandom reader random case analysis\n")
MS_DEN_DIFF_FOM_RRRC <- retMS$msTR+max(retMS$msTC - retMS$msTRC,0) # den of Eqn. (9.23)
F_DBMH <- retMS$msT / MS_DEN_DIFF_FOM_RRRC # Eqn. (9.23)
ndf <- (I-1)
ddfH <- MS_DEN_DIFF_FOM_RRRC^2/(retMS$msTR^2/((I-1)*(J-1))) # Eqn. (9.22)
cat("Hillis ddfH = ", ddfH, "\n")
FCrit <- qf(1 - alpha, ndf, ddfH)
cat("F statistic is ", F_DBMH, ", the critical value of F is ", FCrit, "\n")  
pValueH <- 1 - pf(F_DBMH, ndf, ddfH);cat("p-value is ", pValueH, "\n")

trtMeans <- array(dim = I)
for (i in 1:I) trtMeans[i] <- mean(FOM[i,])
trtDiff <- array(dim = c(I,I))
for (i1 in 1:(I-1)) {    
  for (i2 in (i1+1):I) {
    trtDiff[i1,i2] <- trtMeans[i1]- trtMeans[i2]    
  }
}
trtDiff <- trtDiff[!is.na(trtDiff)]

std_DIFF_FOM_RRRC <- sqrt(2*MS_DEN_DIFF_FOM_RRRC/J/K)
nDiffs <- I*(I-1)/2
CI_DIFF_FOM_RRRC <- array(dim = c(nDiffs, 3))
for (i in 1 : nDiffs) {
  CI_DIFF_FOM_RRRC[i,1] <- trtDiff[i]
  CI_DIFF_FOM_RRRC[i,2] <- qt(alpha/2,df = ddfH)*std_DIFF_FOM_RRRC + trtDiff[i]
  CI_DIFF_FOM_RRRC[i,3] <- qt(1-alpha/2,df = ddfH)*std_DIFF_FOM_RRRC + trtDiff[i]
  cat("mean diff is ", CI_DIFF_FOM_RRRC[i,1], 
      " and 95% CI is ", CI_DIFF_FOM_RRRC[i,2], CI_DIFF_FOM_RRRC[i,3], "\n")
}

cat("\nFixed reader random case analysis\n")
MS_DEN_DIFF_FOM_FRRC <- retMS$msTC
FDbmFR <- retMS$msT / MS_DEN_DIFF_FOM_FRRC
ndf <- (I-1)
ddf <- (I-1)*(K-1)
cat("ddf = ", ddf, "\n")
FCrit <- qf(1 - alpha, ndf, ddf);cat("F statistic is ", FDbmFR, 
                                     "and critical value of F is ", FCrit, "\n")  
pValue <- 1 - pf(FDbmFR, ndf, ddf);cat("p-value is ", pValue, "\n")

std_DIFF_FOM_FRRC <- sqrt(2*MS_DEN_DIFF_FOM_FRRC/J/K)
nDiffs <- I*(I-1)/2
CI_DIFF_FOM_FRRC <- array(dim = c(nDiffs, 3))
for (i in 1 : nDiffs) {
  CI_DIFF_FOM_FRRC[i,1] <- trtDiff[i]
  CI_DIFF_FOM_FRRC[i,2] <- qt(alpha/2,df = ddf)*std_DIFF_FOM_FRRC + trtDiff[i]
  CI_DIFF_FOM_FRRC[i,3] <- qt(1-alpha/2,df = ddf)*std_DIFF_FOM_FRRC + trtDiff[i]
  cat("mean diff is ", CI_DIFF_FOM_FRRC[i,1], 
      " and 95% CI is ", CI_DIFF_FOM_FRRC[i,2], CI_DIFF_FOM_FRRC[i,3], "\n")
}

cat("\nRandom reader fixed case analysis\n")
FDbmFC <- retMS$msT / retMS$msTR
ndf <- (I-1)
ddf <- (I-1)*(J-1)
cat("ddf = ", ddf, "\n")
FCrit <- qf(1 - alpha, ndf, ddf);cat("F statistic is ", FDbmFC, 
                                     "and critical value of F is ", FCrit, "\n")  
pValue <- 1 - pf(FDbmFC, ndf, ddf);cat("p-value is ", pValue, "\n")

MS_DEN_DIFF_FOM_RRFC <- retMS$msTR
std_DIFF_FOM_RRFC <- sqrt(2*MS_DEN_DIFF_FOM_RRFC/J/K)
nDiffs <- I*(I-1)/2
CI_DIFF_FOM_RRFC <- array(dim = c(nDiffs, 3))
for (i in 1 : nDiffs) {
  CI_DIFF_FOM_RRFC[i,1] <- trtDiff[i]
  CI_DIFF_FOM_RRFC[i,2] <- qt(alpha/2,df = ddf)*std_DIFF_FOM_RRFC + trtDiff[i]
  CI_DIFF_FOM_RRFC[i,3] <- qt(1-alpha/2,df = ddf)*std_DIFF_FOM_RRFC + trtDiff[i]
  cat("mean diff is ", CI_DIFF_FOM_RRFC[i,1], 
      " and 95% CI is ", CI_DIFF_FOM_RRFC[i,2], CI_DIFF_FOM_RRFC[i,3], "\n")
}

