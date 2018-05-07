# mainFrocVsAfroc.R
rm(list = ls())
library(RJafroc)
library(ggplot2)
library(abind)

seed <- 1;set.seed(seed)
# following parameters do not change between RAD and CAD observers
nu <- 1;lambda <- 1;K1 <- 500;K2 <- 700 
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1)) # do

#########################  CAD ################################
mu <- 1;zeta1 <- -1 # these are CAD parameters, 13.6 plot a, b
mu <- 1;zeta1 <- -Inf # these are CAD parameters 13.6 plot c, d
#mu <- 1;zeta1 <- -1 # these are CAD parameters, 13.7 plot a, b
#mu <- 1;zeta1 <- -Inf # these are CAD parameters, 13.7 plot c, d
cat("constant parameters:", 
    "\nnu = ", nu, 
    "\nlambda = ", lambda, 
    "\nK1 = ", K1, 
    "\nK2 = ", K2,
    "\n")
cat("CAD parameters: \nmu = ", mu,
    "\nzeta1 = ", zeta1 ,"\n")
frocDataCad <- SimulateFrocDataset(
  mu = mu, lambda = lambda, nu = nu, 
  I = 1, J = 1, K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1
)

#########################  RAD ################################
mu <- 1.5;zeta1 <- 1.5 # these are RAD parameters 13.6 a, b
mu <- 1.5;zeta1 <- -Inf # these are RAD parameters, 13.6 plots c,d
mu <- 2;zeta1 <- 2 # these are RAD parameters, plots 13.7, a, b
mu <- 2;zeta1 <- -Inf # these are RAD parameters, plots 13.7, c, d
#mu <- 1.1;zeta1 <- -1 # these are RAD parameters, 13.8 plot a, b
#mu <- 1.1;zeta1 <- -Inf # these are RAD parameters, 13.8 plot c,d

cat("RAD parameters: \nmu = ", mu, 
    "\nzeta1 = ", zeta1 ,"\n")
set.seed(seed)
frocDataRad <- SimulateFrocDataset(
  mu = mu, lambda = lambda, nu = nu, 
  I = 1, J = 1, K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1
)

numNL1 <- dim(frocDataCad$NL)[4]
numNL2 <- dim(frocDataRad$NL)[4]
numNL <- max(numNL1, numNL2) # the max number of NLs in the combined dataset

if (numNL1 < numNL){ # dataset 1 has smaller number of NLs
  NL <- frocDataCad$NL
  # add more -Inf NLs to make the number of NL in two datasets consistent
  NL <- abind(NL, array(-Inf, dim = c(1, 1, K1 + K2, numNL - numNL1)))
  NL <- abind(NL, frocDataRad$NL, along = 2) # combine the two NLs
}else if (numNL2 < numNL){ # dataset 2 has smaller number of NLs
  NL <- frocDataRad$NL
  NL <- abind(NL, array(-Inf, dim = c(1, 1, K1 + K2, numNL - numNL2)))
  NL <- abind(frocDataCad$NL, NL, along = 2)
}else{ # the number of NLs in the two datasets are same, combine them directly
  NL <- frocDataCad$NL
  NL <- abind(NL, frocDataRad$NL, along = 2)
}

LL <- frocDataCad$LL
LL <- abind(LL, frocDataRad$LL, along = 2) # combine the two LLs

attr(NL, "dimnames") <- NULL; attr(LL, "dimnames") <- NULL

frocDataRaw <- Df2RJafrocDataset(NL, LL, lesionNum = Lk2) 

wAfroc <- UtilFigureOfMerit(frocDataRaw);cat("wAfroc = ", wAfroc, "\n")

froc <- PlotEmpiricalOperatingCharacteristics(frocDataRaw, trts= 1, rdrs = c(1, 2), opChType = "FROC")
combinedPlot <- froc$Plot + scale_color_manual(labels = c("CAD", "RAD"), values = c("black","darkgrey"))  +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(0.8,0.1), legend.direction = "horizontal",
        legend.text = element_text(size = 15, face = "bold"),legend.key.width = unit(1.5, "cm")) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) 
combinedPlot$layers[[1]]$aes_params$size <- 2 # line
#combinedPlot$layers[[2]]$aes_params$size <- 5 # points
print(combinedPlot)

afroc <- PlotEmpiricalOperatingCharacteristics(frocDataRaw, trts= 1, rdrs = c(1, 2), opChType = "AFROC")
combinedPlot <- afroc$Plot + scale_color_manual(labels = c("CAD", "RAD"), values = c("black","darkgrey"))  +
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(0.8,0.1), legend.direction = "horizontal",
        legend.text = element_text(size = 15, face = "bold"),legend.key.width = unit(1.5, "cm")) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 
combinedPlot$layers[[1]]$aes_params$size <- 2 # line
#combinedPlot$layers[[2]]$aes_params$size <- 5 # points
print(combinedPlot)
