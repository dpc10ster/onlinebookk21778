DfBinDataset2 <- function(dataset, desiredNumBins = 7) {
  I <- length(dataset$NL[,1,1,1])
  J <- length(dataset$NL[1,,1,1])
  K1 <- length(dataset$NL[1,1,,1])
  K2 <- length(dataset$LL[1,1,,1])
  K1 <- K1 - K2
  fomOrg <- as.matrix(UtilFigureOfMerit(dataset, fom = "Wilcoxon"), nrow = I, ncol = J)
  print(fomOrg)
  cat("mean, sd = ", mean(fomOrg), sd(fomOrg), "\n")
  numZeta <- desiredNumBins - 1
  maxFomij <- array(-1,dim = c(I,J))
  sSave <- array(dim = c(I,J))
  zetasArr <- array(dim = c(I,J,numZeta))
  datasetB <- dataset
  for (i in 1:I){
    for (j in 1:J){ 
      NL <- dataset$NL[i,j,1:K1,1]
      LL <- dataset$LL[i,j,1:K2,1]
      nLlL <- c(NL,LL)
      # need to remove lowest value, as this gives (1,1) point
      candidateZetas <-  sort(unique(nLlL))[-1]
      el <- length(candidateZetas)
      if (el < numZeta) {
        sample <- combn(candidateZetas, el -1)
      } else {
        # if more than 20 candidates, need to trim
        if (el > 20) {
          byDivisor <- 10
          while (1) {
            by <- as.integer(el/byDivisor)
            candidateZetasTrim <- candidateZetas[seq(from = 1, to = el, by = by)]
            sample <- combn(candidateZetasTrim, numZeta)
            if (length(sample[1,]) > 200) {
              byDivisor <- byDivisor - 1
            } else break
          }
        } else sample <- combn(candidateZetas, numZeta)
      }
      for (s in 1:length(sample[1,])) {
        zetas1 <- sort(sample[,s])
        zetas <- c(-Inf,zetas1,+Inf)
        nLlLB <- cut(nLlL, zetas, labels = FALSE, right = FALSE)
        datasetB$NL[i,j,1:K1,1] <- nLlLB[1:K1]
        datasetB$LL[i,j,1:K2,1] <- nLlLB[(K1+1):(K1+K2)]
        fom <- UtilFigureOfMerit(datasetB, fom = "Wilcoxon")[i,j]
        if (fom > maxFomij[i,j]){
          sSave[i,j] <- s
          maxFomij[i,j] <- fom
          zetasArr[i,j,1:length(zetas1)] <- zetas1
        }
      }
      next
    }
  }
  
  datasetB <- dataset
  for (i in 1:I) {
    for (j in 1:J) { 
      NL <- dataset$NL[i,j,1:K1,1]
      LL <- dataset$LL[i,j,1:K2,1]
      z <- zetasArr[i,j,]
      z <- z[!is.na(z)]
      zetas <- c(-Inf,z,Inf)
      nLlLB <- cut(nLlL, zetas, labels = FALSE, right = FALSE)
      datasetB$NL[i,j,1:K1,1] <- nLlLB[1:K1]
      datasetB$LL[i,j,1:K2,1] <- nLlLB[(K1+1):(K1+K2)]
    }
  }
  fom <- UtilFigureOfMerit(datasetB, fom = "Wilcoxon")
  print(fom)
  cat("mean, sd = ", mean(fom), sd(fom), "\n")
  return(datasetB)
}