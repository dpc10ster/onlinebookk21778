source("MySyncFile.R") # IMPORTANT: ALWAYS chk before analysis
source("TrapezoidalAreaROC.R")
orCovariances <- function (JACK, ratings, truth, aucEst = "RocFit")
{
  # aucEst  <-  "Trap" # Temp*DPC
  I <- length(ratings[,1,1])
  J <- length(ratings[1,,1])
  K <- length(ratings[1,1,])
  
  if (JACK) {
    if (aucEst == "Trap"){
      fomArrayjk <- array(dim = c(I, J, K))
      for (k in 1:length(truth)) {
        ratings_jk <- ratings[ , , -k]
        truth_jk <- truth[ -k ]
        fomArrayjk[ , , k] <- TrapezoidalAreaROC(ratings_jk, truth_jk) 
      }
    }else{
      fomArrayjk <- array(dim = c(I, J, K))
      for (i in 1:I){
        for (j in 1:J){
          for (k in 1:length(truth)) {
            ratings_jk <- ratings[ i, j, -k]
            truth_jk <- truth[ -k ]
            K1b <- table(cut(ratings_jk[!as.logical(truth_jk)], c(1:(DesiredNumBins + 1)), right = FALSE))
            K2b <- table(cut(ratings_jk[as.logical(truth_jk)], c(1:(DesiredNumBins + 1)), right = FALSE))
            Kb <- rbind(K1b, K2b)
            tmp <- RocfitR(Kb)
            if (length(tmp) != 1){
              fomArrayjk[ i, j, k] <- tmp$Az
            }
          }
          if (all(is.na(fomArrayjk[ i, j, ]))){
            fomArrayjkTmp <- rep(NA, length(truth))
            for (k in 1:length(truth)) {
              ratings_jk <- ratings[ i, j, -k]
              truth_jk <- truth[ -k ]
              dim(ratings_jk) <- c(1, 1, length(ratings_jk))
              fomArrayjkTmp[k] <- as.numeric(TrapezoidalAreaROC(ratings_jk, truth_jk))
            }
            dummyStop <- 1
            fomArrayjk[i, j, ] <- fomArrayjkTmp
          }else if (anyNA(fomArrayjk[ i, j, ])){
            kNA <- which(is.na(fomArrayjk[ i, j, ]))
            fomArrayjk[ i, j, kNA ] <- mean(fomArrayjk[ i, j, -kNA ])
          }
        }
      }
    }
    cov <- Cov1Cov2Cov3(fomArrayjk)
    ##############################
    var <- cov$Var * (K - 1)^2/K
    cov1 <- cov$Cov1 * (K - 1)^2/K
    cov2 <- cov$Cov2 * (K - 1)^2/K
    cov3 <- cov$Cov3 * (K - 1)^2/K
    # We missed these four lines before, 
    # but new simulator still worked.
    # Very strange.
    ##############################
    cov2Single <- cov$cov2Single
    varSingle <- cov$varSingle
  } else {
    stop(" BS not implemented yet")
  }
  
  varEchRdr <- vector(length = J)
  cov1EchRdr <- vector(length = J)
  fomArraySinglejk <- array(dim = c(I, 1, K)) 
  for (j in 1:J) {
    fomArraySinglejk[,1,] <- fomArrayjk[, j,]
    Cov <- Cov1Cov2Cov3(fomArraySinglejk)
    varEchRdr[j] = Cov$Var*(K-1)^2/K # see paper by Efron and Stein
    cov1EchRdr[j] = Cov$Cov1*(K-1)^2/K
  }
  
  if (aucEst == "Trap"){
    fomArray <- TrapezoidalAreaROC (ratings, truth) 
  }else{
    fomArray <- apply(fomArrayjk, c(1, 2), mean)
    for (i in 1:I){
      for (j in 1:J){
        K1b <- table(cut(ratings[i, j, ][!as.logical(truth)], c(1:(DesiredNumBins + 1)), right = FALSE))
        K2b <- table(cut(ratings[i, j, ][as.logical(truth)], c(1:(DesiredNumBins + 1)), right = FALSE))
        Kb <- rbind(K1b, K2b)
        fitAuc <- RocfitR(Kb)
        if (length(fitAuc) != 1){
          fomArray[i, j] <- fitAuc$Az
        }
      }
    }
  }
  
  f <- NA
  ddf <- NA
  pValue <- NA
  if (nhCondition == TRUE){
    fomMean <- mean(fomArray)
    msT <- 0
    for (i in 1:I) {
      msT <- msT + (mean(fomArray[i, ]) - fomMean)^2
    }
    msT <- J * msT/(I - 1)
    
    msTR <- 0
    for (i in 1:I) {
      for (j in 1:J) {
        msTR <- msTR + (fomArray[i, j] - mean(fomArray[i, ]) - mean(fomArray[, j]) + fomMean)^2
      }
    }
    msTR <- msTR/((J - 1) * (I - 1))
    
    msNum <- msT
    
    msDenRRRC <- msTR + max(J * (cov2 - cov3), 0)
    f <- msNum/msDenRRRC
    ddf <- msDenRRRC^2/(msTR^2/((I - 1) * (J - 1)))
    pValue <- 1 - pf(f, I - 1, ddf)
  }
  
  return (list(
    fomArray = fomArray,
    var = var,
    cov1 = cov1,
    cov2 = cov2,
    cov3 = cov3,
    cov2Single = cov2Single,
    varSingle = varSingle,
    varEchRdr = varEchRdr,
    cov1EchRdr = cov1EchRdr,
    f = f,
    ddf = ddf,
    pValue = pValue
  ))
}

Cov1Cov2Cov3 <- function (PseudovalueMatrix) 
{ 
  I <- dim(PseudovalueMatrix)[1]
  J <- dim(PseudovalueMatrix)[2]
  Covariance <- array(dim = c(I, I, J, J))
  
  for (i in 1:I){
    for (ip in 1:I){
      for ( j in 1:J) {   
        for ( jp in 1:J) {           
          Covariance[i, ip, j, jp] <- cov(PseudovalueMatrix[i, j, ], PseudovalueMatrix[ip, jp, ])          
        }
      }
    }
  }  
  
  Var <- 0
  count <- 0
  for (i in 1:I){    
    for (j in 1:J) {      
      Var <- Var + Covariance[i, i, j, j] 
      count <- count + 1
    }
  }
  Var <- Var / count 
  
  varSingle <- array(0, dim = I)
  for (i in 1:I){    
    count <- 0
    for (j in 1:J) {      
      varSingle[i] <- varSingle[i] + Covariance[i, i, j, j] 
      count <- count + 1
    }
    varSingle[i] <- varSingle[i] / count   
  }
  
  Cov1 <- 0
  count <- 0
  for (i in 1:I){    
    for (ip in 1:I){
      for (j in 1:J) {      
        if (ip != i){
          Cov1 <- Cov1 + Covariance[i, ip, j, j] 
          count <- count + 1
        }
      }
    }
  }  
  Cov1 <- Cov1 / count 
  
  Cov2 <- 0
  count <- 0
  for (i in 1:I){    
    for (j in 1:J) {      
      for (jp in 1:J){
        if (j != jp){
          Cov2 <- Cov2 + Covariance[i, i, j, jp] 
          count <- count + 1
        }
      }
    }
  }  
  #Cov2 <- Cov2 / (I*J*(J-1)) # OK, DPC
  Cov2 <- Cov2 / count 
  
  # this is needed for single modality stats
  cov2Single <- array(0, dim = I)
  for (i in 1:I){    
    count <- 0
    for (j in 1:J) {      
      for (jp in 1:J){
        if (j != jp){
          cov2Single[i] <- cov2Single[i] + Covariance[i, i, j, jp] 
          count <- count + 1
        }
      }
    }
    cov2Single[i] <- cov2Single[i] / count 
  }  
  
  Cov3 <- 0
  count <- 0
  for (i in 1:I){
    for (ip in 1:I){
      if (i != ip){
        for (j in 1:J) {
          for (jp in 1:J){
            if (j != jp) {
              Cov3 <- Cov3 + Covariance[i, ip, j, jp] 
              count <- count + 1
            }
          }
        }
      }
    }
  }
  
  #Cov3 <- Cov3 / (I*(I-1)*J*(J-1)) # not OK; general advice; better to let computer do the thinking
  Cov3 <- Cov3 / count
  
  return (list (    
    Var = Var,
    Cov1 = Cov1,
    Cov2 = Cov2,
    Cov3 = Cov3,
    varSingle = varSingle,
    cov2Single = cov2Single
  ))  
}

