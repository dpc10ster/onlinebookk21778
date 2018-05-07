rm(list = ls())
library(RJafroc)
library(doParallel)
library(doRNG)
library(foreach)
library(mvtnorm)
library(lmf)

stop("fix or delete me")

fileName <- "NICO"
rocData <- loadDataFile(fileName, "./")
retOrh <- UtilVarianceComponents(rocData, "Wilcoxon", method = "ORH")

I <- length(rocData$modalityID)
J <- length(rocData$readerID)

numCmb <- (I + I * (I - 1)/2) * (J + J * (J - 1)/2)
cmbIJIndx <- array(dim = c(numCmb, 4))
count <- 1
for (i1 in 1:I){
  for (i2 in i1:I){
    for (j1 in 1:J){
      for (j2 in j1:J){
        cmbIJIndx[count, ] <- c(i1, i2, j1, j2)
        count <- count + 1
      }
    }
  }
}

retFile <- paste0("RETCORCBM/", fileName)
if (!file.exists(retFile)){
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  retCorCbm <- foreach(s = 1:numCmb, .errorhandling = "pass", .packages = "RJafroc") %dopar% {
    i1 <- cmbIJIndx[s, 1]; i2 <- cmbIJIndx[s, 2]; j1 <- cmbIJIndx[s, 3]; j2 <- cmbIJIndx[s, 4]
    if (i1 == i2 && j1 == j2){
      retTemp <- FitCbmRoc(rocData, i1, j1)
      list(params = c(retTemp$mu, retTemp$mu, retTemp$alpha, retTemp$alpha, 1, 1),
           covMat = retTemp$covMat[[1]])
    }else{
      retTemp <- FitCorCbmRoc(rocData, c(i1, i2), c(j1, j2))
      list(params = c(retTemp$muX, retTemp$muY, retTemp$alphaX, retTemp$alphaY, retTemp$rhoNor, retTemp$rhoAbn2),
           covMat = retTemp$covMat)
    }
  }
  stopCluster(cl)
  save(retCorCbm, file = retFile)
}else{
  load(retFile)
}

params <- unlist(lapply(retCorCbm, "[[", 1))
dim(params) <- c(6, numCmb)
params <- t(params)

params1 <- rep(0, 2)
countParams1 <- rep(0, 2)
varVector1 <- rep(0, 2)
countVar1 <- rep(0, 2)
corMatr1 <- array(0, dim = c(2, 2))
countMatr1 <- array(0, dim = c(2, 2))

params2 <- rep(0, 6)
countParams2 <- rep(0, 6)
varVector2 <- rep(0, 6)
countVar2 <- rep(0, 6)
corMatr2 <- array(0, dim = c(6, 6))
countMatr2 <- array(0, dim = c(6, 6))

for (s in 1:numCmb){
  i1 <- cmbIJIndx[s, 1]; i2 <- cmbIJIndx[s, 2]; j1 <- cmbIJIndx[s, 3]; j2 <- cmbIJIndx[s, 4]
  if (i1 == i2 && j1 == j2){
    params1 <- params1 + retCorCbm[[s]]$params[c(1, 3)]
    countParams1 <- countParams1 + 1
    covMat <- retCorCbm[[s]]$covMat
    varTemp <- diag(covMat)
    notIsNA <- !is.na(varTemp)
    countVar1 <- countVar1 + notIsNA
    varTemp[is.na(varTemp)] <- 0
    varVector1 <- varVector1 + varTemp
    
    corMatrTemp <- cov2cor(covMat)
    notIsNA <- !is.na(corMatrTemp)
    countMatr1 <- countMatr1 + notIsNA
    corMatrTemp[is.na(corMatrTemp)] <- 0
    corMatr1 <- corMatr1 + corMatrTemp
  }else{
    params2 <- params2 + retCorCbm[[s]]$params
    countParams2 <- countParams2 + 1
    covMat <- retCorCbm[[s]]$covMat
    varTemp <- diag(covMat)
    notIsNA <- !is.na(varTemp)
    countVar2 <- countVar2 + notIsNA
    varTemp[is.na(varTemp)] <- 0
    varVector2 <- varVector2 + varTemp
    
    corMatrTemp <- cov2cor(covMat)
    notIsNA <- !is.na(corMatrTemp)
    countMatr2 <- countMatr2 + notIsNA
    corMatrTemp[is.na(corMatrTemp)] <- 0
    corMatr2 <- corMatr2 + corMatrTemp
    
    temp <- retCorCbm[[s]]$params
    temp[1:4] <- temp[c(2, 1, 4, 3)]
    params2 <- params2 + temp
    countParams2 <- countParams2 + 1
    
    covMat <- retCorCbm[[s]]$covMat
    covMat[1:4, ] <- covMat[c(2, 1, 4, 3), ]
    covMat[ , 1:4] <- covMat[ , c(2, 1, 4, 3)]
    varTemp <- diag(covMat)
    notIsNA <- !is.na(varTemp)
    countVar2 <- countVar2 + notIsNA
    varTemp[is.na(varTemp)] <- 0
    varVector2 <- varVector2 + varTemp
    
    corMatrTemp <- cov2cor(covMat)
    notIsNA <- !is.na(corMatrTemp)
    countMatr2 <- countMatr2 + notIsNA
    corMatrTemp[is.na(corMatrTemp)] <- 0
    corMatr2 <- corMatr2 + corMatrTemp
  }
}
params1 <- params1/countParams1
varVector1 <- varVector1/countVar1
corMatr1 <- corMatr1/countMatr1
covMatr1 <- array(dim = c(2, 2))
diag(covMatr1) <- varVector1
covMatr1[1, 2] <- covMatr1[2, 1] <- sqrt(varVector1[1] * varVector1[2]) * corMatr1[1, 2]

params2 <- params2/countParams2
params2[1:2] <- (params2[1:2] + params1[1])/2
params2[3:4] <- (params2[3:4] + params1[2])/2
varVector2 <- varVector2/countVar2
corMatr2 <- corMatr2/countMatr2
corMatr2[1, 3] <- corMatr2[3, 1] <- (corMatr2[1, 3] + corMatr1[1, 2])/2
corMatr2[2, 4] <- corMatr2[4, 2] <- (corMatr2[2, 4] + corMatr1[1, 2])/2
covMatr2 <- array(dim = c(6, 6))
diag(covMatr2) <- varVector2
varVector2[1:2] <- (varVector2[1:2] + varVector1[1])/2
varVector2[3:4] <- (varVector2[3:4] + varVector1[2])/2
for (iRow in 1:6){
  for (iCol in iRow:6){
    covMatr2[iRow, iCol] <- covMatr2[iCol, iRow] <- sqrt(varVector2[iRow] * varVector2[iCol]) * corMatr2[iRow, iCol]
  }
}


simuJ <- 10; K1 <- 120; K2 <- 80
seed <- 1
cl <- makeCluster(detectCores())
registerDoParallel(cl)
S <- 2000
retSimu <- foreach(s = 1:S, .combine = "rbind", .packages = c("mvtnorm", "RJafroc", "lmf"), .options.RNG = seed) %dorng%{
  # muAlpha <- rmvnorm(simuJ, mean = params1[1:2], sigma = covMatr1)
  simuMu <- rep(0, simuJ)
  simuAlpha <- rep(0, simuJ)
  rhoNorMatr <- array(dim = c(simuJ, simuJ))
  rhoAbn2Matr <- array(dim = c(simuJ, simuJ))
  countJ <- rep(0, simuJ)
  for (j1 in 1:simuJ){
    for (j2 in j1:simuJ){
      if (j1 == j2){
        rhoNorMatr[j1, j2] <- 1
        rhoAbn2Matr[j1, j2] <- 1
      }else{
        paramsTemp <- rmvnorm(1, mean = params2, sigma = covMatr2)
        rhoNorMatr[j1, j2] <- rhoNorMatr[j2, j1] <- paramsTemp[5]
        rhoAbn2Matr[j1, j2] <- rhoAbn2Matr[j2, j1] <- paramsTemp[6]
        simuMu[c(j1, j2)] <- simuMu[c(j1, j2)] + paramsTemp[c(1, 2)]
        simuAlpha[c(j1, j2)] <- simuAlpha[c(j1, j2)] + paramsTemp[c(3, 4)]
        countJ[c(j1, j2)] <- countJ[c(j1, j2)] + 1
      }
    }
  }
  simuMu <- simuMu/countJ
  simuAlpha <- simuAlpha/countJ
  simuAlpha[simuAlpha > 1] <- 1
  conditionArray <- rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1))
  for (j in 3:simuJ){
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
    sigmaTemp <- array(dim = c(simuJ, simuJ))
    for (iRow in 1:simuJ){
      for (iCol in iRow:simuJ){
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
  dim(z1) <- c(1, simuJ, K1, 1)
  dim(z2) <- c(1, simuJ, K2, 1)
  stop("check me")
  simuData <- Df2RJafrocDataset(z1, z2)
  retOrhSimu <- UtilVarianceComponents(simuData, fom = "Wilcoxon", method = "ORH")
  ret <- StSignificanceTestingCadVsRadiologists(simuData, fom = "Wilcoxon")
  c(retOrhSimu$cov2, retOrhSimu$varEps, ret$pVal)
}
stopCluster(cl)
options(digits = 6)
means <- colMeans(retSimu)
cov2Mean <- means[1]
varMean <- means[2]

cat(simuJ, ",", K1, ",", K2, "\n")
cat("cov2 and var of Orig. data are:", retOrh$cov2, retOrh$varEps, "\n")
cat("cov2 and var of Simu. data are:", cov2Mean, varMean, "\n")
cat("95% CI of cov2 is: (", quantile(retSimu[ , 1], 0.025), ",", quantile(retSimu[ , 1], 0.975), ")\n")
cat("95% CI of var is: (", quantile(retSimu[ , 2], 0.025), ",", quantile(retSimu[ , 2], 0.975), ")\n")
cat("Rej. rate is", sum(retSimu[ , 3] < 0.05)/S, "\n\n")



