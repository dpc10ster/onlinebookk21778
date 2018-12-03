rm(list = ls()) # mainRocFitR.R 
library("numDeriv")
library(RJafroc)

source("Transforms.R")
source("LL.R")
source("RocOperatingPointsFromRatingsTable.R")
source("VarianceAz.R")
source("ChisqrGoodnessOfFit.R");source("UtilGoodnessOfFit.R")
source("CombineBins.R")

options(digits = 4)
# clamps on range of allowed values
a_min <<- 0.001;a_max <<- 6
b_min <<- 0.001;b_max <<- 6

K1 <- c(30,19,8,2,1) # this is the observed data!
#K1 <- c(30,19,8,7,5) # this is the cheated data!
K2 <- c(5,6,5,12,22) # this is the observed data!

# initial estimates of a and b parameters
ret <- RocOperatingPointsFromRatingsTable (K1, K2)
FPF <- ret$FPF; TPF <- ret$TPF

phiInvFpf <- qnorm(FPF)
phiInvTpf <- qnorm(TPF)
# straight line fit method of estimating a and b
fit <- lm(phiInvTpf~phiInvFpf)
# these is the initial estimate of a and b
a <- fit$coefficients[[1]]
b <- fit$coefficients[[2]]

# thresholds can be estimated by 
# applying inverse function to Eqn. xx and 
# solving for zeta
zetaIniFpf <- -phiInvFpf
zetaIniTpf <- (a - phiInvTpf)/b
zetaIni <- (zetaIniFpf + zetaIniTpf)/2 # average the two estimates
# apply reverse order to correct the ordering of the cutoffs
zetaIni <- rev(zetaIni)
# to test stability of alg. to guess choice
#zetaIniGuess <- seq(-b, a + 1, length.out = length(K1)-1)

paramIni <- c(a, b, zetaIni) 
# to test stability of alg. to other choices
#paramIni <- c(1, 1, zetaIniGuess)

# use this method to test variation of -LL with parameters
paramIniPrime <- ThetaPrime(paramIni)
LLvalIni <- LL(paramIniPrime, K1, K2)

# this does the actual minimization of -LL
retNlm <- nlm(LL, 
              paramIniPrime, 
              K1 = K1, 
              K2 = K2, 
              stepmax = 0.1)
paramFinal <- Theta(retNlm$estimate)

hess <- hessian(
  LL_theta, 
  paramFinal, 
  method="Richardson", 
  K1 = K1, K2 = K2)
Cov <- solve(hess)
a <- paramFinal[1]
b <- paramFinal[2]
Az <- pnorm(a/sqrt(1+b^2))
StdAz <- sqrt(VarianceAz (a, b,Cov))

cat("initial parameters = \n", paramIni, "\n")
cat("final parameters = \n", paramFinal, "\n")
cat("-LL values, initial", LLvalIni, "\n")
cat("-LL values, final", retNlm$minimum,"\n")
cat("covariance matrix = \n")
for (i in 1:6){
  for (j in 1:6){
    x <- sprintf("%7.3f", Cov[i,j]); cat(x)
  }
  cat("\n")
}
cat("\nAz = ", Az, "StdAz = ", StdAz, "\n")

options(digits = 3)
retChisqInitial <- ChisqrGoodnessOfFit(paramIni,K1,K2)
retChisqFinal <- ChisqrGoodnessOfFit(paramFinal,K1,K2)
if (!anyNA(retChisqInitial)) 
  cat("retChisqInitial p-val = ", retChisqInitial$pVal,"\n")
if (!anyNA(retChisqFinal)) 
  cat("Chisq = ", retChisqFinal$chisq, 
      "\nChisq df = ", retChisqFinal$df, 
      "\nChisq p-val = ", retChisqFinal$pVal, "\n")

# ML estimates according to Eng program
# FINAL VALUES OF PARAMETERS:
#   Procedure converges after 5 iterations.
# A = 1.3204
# B = 0.6075
# Z(K):  0.0077  0.8963  1.5157  2.3967
# LOGL = -141.4354
# 
# VARIANCE-COVARIANCE MATRIX:
# A     0.0656  0.0259  0.0150  0.0128  0.0070 -0.0135
# B     0.0259  0.0254  0.0053 -0.0022 -0.0141 -0.0458
# Z(1)  0.0150  0.0053  0.0260  0.0153  0.0109  0.0042
# Z(2)  0.0128 -0.0022  0.0153  0.0317  0.0276  0.0285
# Z(3)  0.0070 -0.0141  0.0109  0.0276  0.0539  0.0660
# Z(4) -0.0135 -0.0458  0.0042  0.0285  0.0660  0.1664
# 
# SUMMARY OF ROC CURVE:
#   Area = 0.8705
# Std. Dev. (Area) = 0.0378

# output of this program
# initial parameters =  1.328148 0.6292443 0.03702537 0.8931309 1.506143 2.239335 
# final parameters =  1.320453 0.607497 0.007675259 0.8962713 1.515645 2.39671 
# -LL values, initial, final 141.6644 141.4354 
# covariance matrix = 
#   [,1]         [,2]        [,3]         [,4]         [,5]        [,6]
# [1,]  0.065222452  0.025075693 0.014901302  0.012449852  0.007256201 -0.01067237
# [2,]  0.025075693  0.024258552 0.005134078 -0.002321735 -0.014328108 -0.04286043
# [3,]  0.014901302  0.005134078 0.025927200  0.015224066  0.010850408  0.00464582
# [4,]  0.012449852 -0.002321735 0.015224066  0.031623192  0.027831803  0.02855593
# [5,]  0.007256201 -0.014328108 0.010850408  0.027831803  0.055993440  0.06697530
# [6,] -0.010672373 -0.042860433 0.004645820  0.028555929  0.066975300  0.15851205
# Az =  0.8695184 StdAz =  0.0376075 


