rm(list = ls()) # mainExamples.R
require(RJafroc)
seed <- 1;set.seed(seed)
K1 <- 4;K2 <- 4
maxLL <- 2;Lk2 <- ceiling(runif(K2) * maxLL) 
mu <- 2;lambda <- 1;nu <- 1 ;zeta1 <- -1
frocData <- SimulateFrocDataset(
  mu = mu, 
  lambda = lambda, 
  nu = nu, 
  I = 1, 
  J = 1,
  K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1)

FP <- apply(frocData$NL, 3, max);FP <- FP[1:K1]
frocData$lesionWeight[3, ] <- c(0.6, 0.4)
frocData$lesionWeight[4, ] <- c(0.4, 0.6)

plotAfroc <- PlotEmpiricalOperatingCharacteristics(
  frocData,
  trts = 1, 
  rdrs = 1, 
  opChType = "AFROC")
print(plotAfroc$Plot)

plotwAfroc <- PlotEmpiricalOperatingCharacteristics(
  frocData,
  trts = 1, 
  rdrs = 1, 
  opChType = "wAFROC")
print(plotwAfroc$Plot)

cat("AFROC AUC = ", UtilFigureOfMerit(frocData, FOM = "AFROC"),"\n")
cat("wAFROC AUC = ", UtilFigureOfMerit(frocData, FOM = "wAFROC"),"\n")