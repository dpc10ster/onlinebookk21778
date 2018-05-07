calibrateMutivariateSimulator <- function(rocData)
{
  retOrh <- StSignificanceTestingCadVsRadiologists (rocData, fom = "Wilcoxon")
  
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
  fileName <- "NICO"
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
  
  paramsCAD <- rep(0, 6)
  countParamsCAD <- rep(0, 6)
  varVectorCAD <- rep(0, 6)
  countVarCAD <- rep(0, 6)
  corMatrCAD <- array(0, dim = c(6, 6))
  countMatrCAD <- array(0, dim = c(6, 6))
  
  for (s in 1:numCmb){
    i1 <- cmbIJIndx[s, 1]; i2 <- cmbIJIndx[s, 2]; j1 <- cmbIJIndx[s, 3]; j2 <- cmbIJIndx[s, 4]
    if (i1 == i2 && j1 == j2){
      if (j1 == 1){
        paramsCAD[1] <- paramsCAD[1] + retCorCbm[[s]]$params[1]
        paramsCAD[3] <- paramsCAD[3] + retCorCbm[[s]]$params[3]
        countParamsCAD[c(1, 3)] <- countParamsCAD[c(1, 3)] + 1
        
        covMat <- retCorCbm[[s]]$covMat
        varTemp <- diag(covMat)
        notIsNA <- !is.na(varTemp)
        countVarCAD[c(1, 3)] <- countVarCAD[c(1, 3)] + notIsNA
        varTemp[is.na(varTemp)] <- 0
        varVectorCAD[c(1, 3)] <- varVectorCAD[c(1, 3)] + varTemp
        assign("last.warning", NULL, envir = baseenv())
        corMatrTemp <- cov2cor(covMat)
        notIsNA <- !is.na(corMatrTemp)
        countMatrCAD[1, 1] <- countMatrCAD[1, 1] + notIsNA[1, 1]
        countMatrCAD[3, 3] <- countMatrCAD[3, 3] + notIsNA[2, 2]
        countMatrCAD[1, 3] <- countMatrCAD[1, 3] + notIsNA[1, 2]
        countMatrCAD[3, 1] <- countMatrCAD[3, 1] + notIsNA[2, 1]
        corMatrTemp[is.na(corMatrTemp)] <- 0
        corMatrCAD[1, 1] <- corMatrCAD[1, 1] + corMatrTemp[1, 1]
        corMatrCAD[3, 3] <- corMatrCAD[3, 3] + corMatrTemp[2, 2]
        corMatrCAD[1, 3] <- corMatrCAD[1, 3] + corMatrTemp[1, 2]
        corMatrCAD[3, 1] <- corMatrCAD[3, 1] + corMatrTemp[2, 1]
      }else{
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
      }
    }else{
      if (j1 == 1){
        paramsCAD <- paramsCAD + retCorCbm[[s]]$params
        countParamsCAD <- countParamsCAD + 1
        covMat <- retCorCbm[[s]]$covMat
        varTemp <- diag(covMat)
        notIsNA <- !is.na(varTemp)
        countVarCAD <- countVarCAD + notIsNA
        varTemp[is.na(varTemp)] <- 0
        varVectorCAD <- varVectorCAD + varTemp
        
        corMatrTemp <- cov2cor(covMat)
        notIsNA <- !is.na(corMatrTemp)
        countMatrCAD <- countMatrCAD + notIsNA
        corMatrTemp[is.na(corMatrTemp)] <- 0
        corMatrCAD <- corMatrCAD + corMatrTemp
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
  }
  
  paramsCAD <- paramsCAD/countParamsCAD
  varVectorCAD <- varVectorCAD/countVarCAD
  corMatrCAD <- corMatrCAD/countMatrCAD
  covMatrCAD <- array(dim = c(6, 6))
  diag(covMatrCAD) <- varVectorCAD
  for (iRow in 1:6){
    for (iCol in iRow:6){
      covMatrCAD[iRow, iCol] <- covMatrCAD[iCol, iRow] <- sqrt(varVectorCAD[iRow] * varVectorCAD[iCol]) * corMatrCAD[iRow, iCol]
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
  
  return(list = list(
    retOrh = retOrh,
    paramsCAD = paramsCAD,
    covMatrCAD = covMatrCAD,
    params2 = params2,
    covMatr2 = covMatr2
  ))
  
}


cov2cor <- function(V)
  {
  assign("last.warning", NULL, envir = baseenv())
  p <- (d <- dim(V))[1L]
    if (!is.numeric(V) || length(d) != 2L || p != d[2L]) 
      stop("'V' is not a square numeric matrix")
    Is <- sqrt(1/diag(V))
    if (any(!is.finite(Is))) {
      cat("diag(.) had 0 or NA entries; non-finite result is doubtful")
      cat("I am stuck here")
      }
    r <- V
    r[] <- Is * V * rep(Is, each = p)
    r[cbind(1L:p, 1L:p)] <- 1
    r
  }
