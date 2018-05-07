# mainCadVsRadCalibValidate.R
rm(list = ls())
library(RJafroc)
library(doParallel)
library(doRNG)
library(foreach)
library(mvtnorm)
library(lmf)
source("calibrateMutivariateSimulator.R")
options(digits = 3)
f <- 9 #NICO
rocData <- get(sprintf("dataset%02d", f)) # the datasets already exist as R objects
retOrgAnalysis <- StSignificanceTestingCadVsRadiologists (rocData, fom = "Wilcoxon")
I <- length(rocData$NL[,1,1,1]);J <- length(rocData$NL[1,,1,1]);
retCalibSimParms <- calibrateMutivariateSimulator(rocData)

retOrh <-  retCalibSimParms$retOrh
paramsCAD = retCalibSimParms$paramsCAD
covMatrCAD = retCalibSimParms$covMatrCAD
params2 = retCalibSimParms$params2
covMatr2 = retCalibSimParms$covMatr2

# NH 
paramsCAD[c(1, 3)] <- (paramsCAD[c(2, 4)] + params2[c(1, 3)])/2
simuJArr <-  c(4, 10, 5, 5, 10, 10, 10, 10, 5,  15)
simuK1Arr <- c(120, 80, 25,50,25, 50, 80, 100,200,200)
simuK2Arr <- c(80, 120,25,50,25, 50, 120,100,200,200)

binning <- TRUE
for (q in 1:length(simuK1Arr)){
  J <- simuJArr[q]; K1 <- simuK1Arr[q]; K2 <- simuK2Arr[q];seed <- NULL; seed <- 1;set.seed(seed)
  cat("seed =", seed, ", binning=", binning, ", J=", J, ", K1=",K1,", K2=",K2,"\n")
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  S <- 2000
  retSimu <- foreach(s = 1:S, .combine = "rbind", .packages = c("mvtnorm", "RJafroc", "lmf"), .options.RNG = seed) %dorng%{
    simuMu <- rep(0, J)
    simuAlpha <- rep(0, J)
    rhoNorMatr <- array(0, dim = c(J, J))
    rhoAbn2Matr <- array(0, dim = c(J, J))
    countJ <- rep(0, J)
    for (j1 in 1:J){
      for (j2 in 1:J){
        if (j1 == j2){
          if (j1 == 1){
            covMatTemp <- rbind(c(covMatrCAD[1, 1], covMatrCAD[1, 3]), c(covMatrCAD[3, 1], covMatrCAD[3, 3]))
            paramsMeanTemp <- paramsCAD[c(1, 3)]
          }else{
            covMatTemp <- rbind(c(covMatr2[1, 1], covMatr2[1, 3]), c(covMatr2[3, 1], covMatr2[3, 3]))
            paramsMeanTemp <- params2[c(1, 3)]
          }
          paramsTemp <- rmvnorm(1, mean = paramsMeanTemp, sigma = nearPD(covMatTemp))
          simuMu[j1] <- simuMu[j1] + paramsTemp[1]
          simuAlpha[j1] <- simuAlpha[j1] + paramsTemp[2]
          countJ[j1] <- countJ[j1] + 1
          rhoNorMatr[j1, j2] <- 2
          rhoAbn2Matr[j1, j2] <- 2
        }else{
          if (j1 == 1){
            paramsTemp <- rmvnorm(1, mean = paramsCAD, sigma = covMatrCAD)
            rhoNorMatr[j1, j2] <- rhoNorMatr[j2, j1] <- rhoNorMatr[j2, j1] + paramsTemp[5]
            rhoAbn2Matr[j1, j2] <- rhoAbn2Matr[j2, j1] <- rhoAbn2Matr[j2, j1] + paramsTemp[6]
            simuMu[c(j1, j2)] <- simuMu[c(j1, j2)] + paramsTemp[c(1, 2)]
            simuAlpha[c(j1, j2)] <- simuAlpha[c(j1, j2)] + paramsTemp[c(3, 4)]
            countJ[c(j1, j2)] <- countJ[c(j1, j2)] + 1
          }else if (j2 == 1){
            paramsCADTemp <- paramsCAD[c(2, 1, 4, 3, 5, 6)]
            covMatrCADTemp <- covMatrCAD[c(2, 1, 4, 3, 5, 6), ]
            covMatrCADTemp <- covMatrCADTemp[ , c(2, 1, 4, 3, 5, 6)]
            paramsTemp <- rmvnorm(1, mean = paramsCADTemp, sigma = covMatrCADTemp)
            rhoNorMatr[j1, j2] <- rhoNorMatr[j2, j1] <- rhoNorMatr[j2, j1] + paramsTemp[5]
            rhoAbn2Matr[j1, j2] <- rhoAbn2Matr[j2, j1] <- rhoAbn2Matr[j2, j1] + paramsTemp[6]
            simuMu[c(j1, j2)] <- simuMu[c(j1, j2)] + paramsTemp[c(1, 2)]
            simuAlpha[c(j1, j2)] <- simuAlpha[c(j1, j2)] + paramsTemp[c(3, 4)]
            countJ[c(j1, j2)] <- countJ[c(j1, j2)] + 1
          }else{
            paramsTemp <- rmvnorm(1, mean = params2, sigma = covMatr2)
            rhoNorMatr[j1, j2] <- rhoNorMatr[j2, j1] <- rhoNorMatr[j2, j1] + paramsTemp[5]
            rhoAbn2Matr[j1, j2] <- rhoAbn2Matr[j2, j1] <- rhoAbn2Matr[j2, j1] + paramsTemp[6]
            simuMu[c(j1, j2)] <- simuMu[c(j1, j2)] + paramsTemp[c(1, 2)]
            simuAlpha[c(j1, j2)] <- simuAlpha[c(j1, j2)] + paramsTemp[c(3, 4)]
            countJ[c(j1, j2)] <- countJ[c(j1, j2)] + 1
          }
        }
      }
    }
    simuMu <- simuMu/countJ
    simuAlpha <- simuAlpha/countJ
    simuAlpha[simuAlpha > 1] <- 1
    rhoNorMatr <- rhoNorMatr/2
    rhoAbn2Matr <- rhoAbn2Matr/2
    rhoNorMatr[rhoNorMatr > 1] <- 1 # samples can exceed one
    rhoAbn2Matr[rhoAbn2Matr > 1] <- 1
 
    conditionArray <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
    for (j in 3:J){
      conditionArray <- rbind(cbind(conditionArray, 0), cbind(conditionArray, 1))
    }
    nCond <- nrow(conditionArray)
    probVector <- rep(NA, nCond)
    for (l in 1:nCond){
      visibleProb <- simuAlpha[which(conditionArray[l, ] == 1)]
      invisibleProb <- 1 - simuAlpha[which(conditionArray[l, ] == 0)]
      probVector[l] <- prod(c(visibleProb, invisibleProb)) 
    }
    
    z1 <- t(rmvnorm(n = K1, sigma = nearPD(rhoNorMatr)))
    K2Sample <- sample(nCond, size = K2, replace = TRUE, prob = probVector)
    uniqueCond <- unique(K2Sample)
    z2 <- NULL
    for (l in uniqueCond){
      k <- sum(K2Sample == l)
      muTemp <- conditionArray[l, ]
      muTemp[muTemp == 1] <- simuMu[muTemp == 1]
      sigmaTemp <- array(dim = c(J, J))
      for (iRow in 1:J){
        for (iCol in iRow:J){
          if (iRow == iCol){
            sigmaTemp[iRow, iCol] <- 1
          }else{
            if (sum(conditionArray[l, iRow] == 1, conditionArray[l, iCol] == 1) == 0){
              sigmaTemp[iRow, iCol] <- sigmaTemp[iCol, iRow] <- rhoNorMatr[iCol, iRow]
            }else if (sum(conditionArray[l, iRow] == 1, conditionArray[l, iCol] == 1) == 1){
              sigmaTemp[iRow, iCol] <- sigmaTemp[iCol, iRow] <- (rhoNorMatr[iCol, iRow] + rhoAbn2Matr[iCol, iRow])/2
            }else if (sum(conditionArray[l, iRow] == 1, conditionArray[l, iCol] == 1) == 2){
              sigmaTemp[iRow, iCol] <- sigmaTemp[iCol, iRow] <- rhoAbn2Matr[iCol, iRow]
            }
          }
        }
      }
      z2Temp <- t(rmvnorm(n = k, mean = muTemp, sigma = nearPD(sigmaTemp)))
      z2 <- cbind(z2, z2Temp)
    }
    dim(z1) <- c(J, K1)
    dim(z2) <- c(J, K2)
    simuData <- Df2RJafrocDataset(z1, z2)
    if (binning){
      simuData <- DfBinDataset(simuData, desiredNumBins = 5)
    }
    ret <- StSignificanceTestingCadVsRadiologists(simuData, fom = "Wilcoxon")
    c(ret$cov2*((K1+K2)/200), ret$var*((K1+K2)/200), ret$pVal)# !dpc
  }
  stopCluster(cl)

  means <- colMeans(retSimu)
  cov2Mean <- means[1]
  varMean <- means[2]
  
  cat("cov2 and var of Orig. data are:", retOrh$cov2, retOrh$varEps, "\n")
  cat("cov2 and var of Simu. data are:", cov2Mean, varMean, "\n")
  lower <- quantile(retSimu[ , 1], 0.025)
  upper <- quantile(retSimu[ , 1], 0.975)
  cat("95% CI of cov2 is: (", lower, ",", upper, ")\n")
  if (cov2Mean >= lower && cov2Mean <= upper){
    cat("cov2 included?", "YES\n")
  }else{
    cat("cov2 included?", "NO\n")
  }
  lower <- quantile(retSimu[ , 2], 0.025)
  upper <- quantile(retSimu[ , 2], 0.975)
  cat("95% CI of var is: (", lower, ",", upper, ")\n")
  if (varMean >= lower && varMean <= upper){
    cat("var included?", "YES\n")
  }else{
    cat("var included?", "No\n")
  }
  cat("Rej. rate is", sum(retSimu[ , 3] < 0.05)/S, "\n\n")
  #stop("temp...")
}


# FALSE	10/80/120	0.00013
# (0.00013, 0.00032)	0.00059
# (0.00050, 0.00082)	0.063
# FALSE	5/25/25	(3.6e-5, 0.00045)	(0.00037, 0.0011)	0.0625
# FALSE	5/50/50	(7.3e-5, 0.00037)	(0.00044, 0.00095)	0.0635
# FALSE	10/25/25	(5.7e-5, 0.00046)	(0.00039, 0.0011)	0.0665
# FALSE	10/50/50	(9.9e-5, 0.00037)	(0.00047, 0.00093)	0.063
# FALSE	10/100/100	(0.00013, 0.00033)	(0.00051, 0.00085)	0.063
# FALSE	5/200/200	(0.00013, 0.00028)	(0.00054, 0.00080)	0.089
# FALSE	15/200/200	(0.00015, 0.0003)	(0.00056,0.00078)	0.062
# 1	TRUE	5/25/25	0.00023
# (5e-05, 0.0005)	0.00085
# (0.00057, 0.0012)	0.04
# 2	TRUE	10/25/25	0.00025
# (0,0.00048)	0.00086
# (0.0006, 0.0012)	0.041
# 
# 1	TRUE	5/50/50	0.00021
# (0, 0.00038)	0.00077
# (0.00056, 0.001)	0.046
# TRUE	10/50/50	0.00023
# (0.00011, 0.00039)	0.00078
# (0.00058, 0.001)	0.047
# TRUE	10/80/120	0.00022*
#   (0.00014, 0.00032)	0.00074
# (0.0006, 0.00089)	0.0465
# TRUE	10/80/120	0.00022*
#   (0.00014, 0.00032)	0.00073
# (0.00060, 0.00088)	0.0485
# TRUE	10/100/100	0.00022*
#   (0.00013, 0.00033)	0.00074
# (0.00059, 0.00089)	0.05
# FALSE	10/100/100	0.00022*
#   (0.00012, 0.00033) 	0.00068*
#   (0.00052, 0.00086)	0.058
