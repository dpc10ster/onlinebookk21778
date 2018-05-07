library(abind)

seed <- 1;set.seed(seed)
nu <- 1;lambda <- 1;K1 <- 50;K2 <- 70 # these parameters do not change between RAD and CAD observers
Lmax <- 2;Lk2 <- floor(runif(K2, 1, Lmax + 1)) # do

mu <- 1;zeta1 <- -1 # these are CAD parameters, plot a, b
frocDataCad <- SimulateFrocDataset(
  mu = mu, lambda = lambda, nu = nu, 
  I = 1, J = 1, K1 = K1, K2 = K2, lesionNum = Lk2, zeta1 = zeta1
)
mu <- 1.5;zeta1 <- 1.5 # these are RAD parameters
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
  NL <- abind(NL, array(-Inf, dim = c(1, 1, K1 + K2, numNL - numNL1))) # add more -Inf NLs to make the number of NL in two datasets consistent
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

frocDataRaw <- Df2RJafrocDataset(NL, LL, lesionNum = Lk2) # convert the the combined NLs and LLs to RJafroc dataset
rocData <- DfFroc2Roc(frocDataRaw)
roc <- PlotEmpiricalOperatingCharacteristics(rocData, trts= 1, rdrs = c(1, 2), opChType = "ROC")
combinedPlot <- froc$Plot + 
  #scale_color_manual(labels = c("CAD", "RAD"), values = c("black","darkgrey")) +
  #theme(legend.position="none") +
  theme(legend.title = element_blank(), legend.position = c(1, 0), 
        legend.key.size = unit(1.5, "lines"), legend.text=element_text(size=20, face = "bold"), legend.direction = "horizontal") +   
  theme(axis.title.y = element_text(size = 25,face="bold"),
        axis.title.x = element_text(size = 30,face="bold"),
        legend.position = c(1,0), legend.direction = "horizontal",
        legend.text = element_text(size = 15, face = "bold"),legend.key.size = unit(1.5, "lines")) 
combinedPlot$layers[[1]]$aes_params$size <- 2 # line
combinedPlot$layers[[2]]$aes_params$size <- 5 # points
print(combinedPlot)
