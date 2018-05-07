CBMAuc <- function(dataset, selectedMod, selectedRdr){
  NL <- dataset$NL
  LL <- dataset$LL
  I <- dim(NL)[1]
  J <- dim(LL)[2]
  
  K <- dim(NL)[3]
  K2 <- dim(LL)[3]
  K1 <- K - K2 
  modalityID <- dataset$modalityID
  readerID <- dataset$readerID
  maxNL <- dim(NL)[4] 
  lesionNum <- dataset$lesionNum
  
  I <- length(selectedMod)
  if (all(is.character(selectedMod)))
    selectedMod <- which( modalityID %in% selectedMod )
  NL <- NL[selectedMod, , 1:K1, ]
  LL <- LL[selectedMod, , , ]
  dim(NL) <- c(I, J, K1, maxNL)
  dim(LL) <- c(I, J, K2, max(lesionNum))
  modalityID <- modalityID[selectedMod]
  
  J <- length(selectedRdr)
  if (all(is.character(selectedRdr)))
    selectedRdr <- which( readerID %in% selectedRdr )
  NL <- NL[ , selectedRdr, , ]
  LL <- LL[ , selectedRdr, , ]
  dim(NL) <- c(I, J, K1, maxNL)
  dim(LL) <- c(I, J, K2, max(lesionNum))
  readerID <- readerID[selectedRdr]
  
  fpBinned <- apply(NL, c(1, 2, 3), max)
  tpBinned <- apply(LL, c(1, 2, 3), max)
  
  nLesDistr <- table(lesionNum)
  if (length(nLesDistr) == 1) {
    nLesDistr <- lesionNum[1]
  }else{
    nLesDistr <- t(rbind(as.numeric(unlist(attr(nLesDistr, "dimnames"))), as.vector(nLesDistr)))
  }
  
  mu <- array(dim = c(I, J))
  alpha <- array(dim = c(I, J))
  AUC <- array(dim = c(I, J))
  gdnss <- array(dim = c(I, J))
  
  maxit <- 10000
  for (i in 1:I){
    for (j in 1:J){
      fpBinnedTable <- table(fpBinned)
      tpBinnedTable <- table(tpBinned)
      if (length(fpBinnedTable) != length(tpBinnedTable) || length(fpBinnedTable) < 4 || length(tpBinnedTable) < 4)
        next
      FPF <- cumsum(rev(fpBinnedTable)) / K1
      TPF <- cumsum(rev(tpBinnedTable)) / K2
      if (FPF[length(FPF) - 1] == 0){
        alphaIni <- TPF[length(TPF) - 1]
        return(0.5 * (1 + alphaIni))
      }else{
        alphaIni <- mean(TPF[1:(length(TPF) - 1)])
        zetaIni <- rep(NA, length(TPF) - 1)
        zetaIni[1] <- qnorm(1 - FPF[length(FPF) - 1])
        y1 <- (TPF[length(TPF) - 1] - (1 - alphaIni) * (1 - FPF[length(FPF) - 1])) / alphaIni
        if (y1 >= 0.999) y1 <- 0.999
        muIni <- zetaIni[1] - qnorm(y1, lower.tail = FALSE)
        if (length(zetaIni) > 1){
          for (z in 2:(length(TPF) - 1)){
            if (z == length(TPF) - 1)
              break
            zetaIni[z] <- optimize(zetaIniOpt, c(zetaIni[z - 1], qnorm(0.999) + muIni), muIni = muIni, alphaIni = alphaIni, yi = TPF[length(TPF) - z])$minimum 
            # zetaIni[z] <- uniroot(zetaIniOpt, c(zetaIni[z - 1], qnorm(0.999) + muIni), muIni = muIni, alphaIni = alphaIni, yi = TPF[length(TPF) - z], tol = .Machine$double.eps^0.5)$root
          }
        }
      }
      # zetaIni[length(TPF) - 1] <- uniroot(zetaIniOpt, c(zetaIni[length(TPF) - 2], qnorm(0.999) + muIni), muIni = muIni, alphaIni = alphaIni, yi = TPF[1], tol = .Machine$double.eps^0.5)$root
      zetaIni[length(TPF) - 1] <- optimize(zetaIniOpt, c(zetaIni[length(TPF) - 2], qnorm(0.999) + muIni), muIni = muIni, alphaIni = alphaIni, yi = TPF[1])$minimum 
      namesVector <- c("mu", "alpha")
      for (t in 1:length(zetaIni))
        namesVector <- c(namesVector, paste0("zetaFwd", t))
      parameters <- c(list(muIni, alphaIni), as.list(zetaIni))
      names(parameters) <- namesVector
      CBMNLLNew <- addArguments(CBMNLL, length(zetaIni))
      if (any(c(fpBinnedTable, tpBinnedTable) == 0)){
        dummy <- 1
      }
      ret <- mle2(CBMNLLNew, start = parameters, method = "Nelder-Mead", data = list(fi = fpBinnedTable, ti = tpBinnedTable), control = list(maxit = maxit))
      mu[i, j] <- ret@coef[1]
      alpha[i, j] <- ret@coef[2]
      AUC[i, j] <- integrate(CBMAucInt, 0, 1, mu = mu[i, j], alpha = alpha[i, j])$value
      gdnss[i, j] <- GdnsFtCBM(ret@coef, fpBinnedTable, tpBinnedTable, K1, K2)
    }
  }
  return(list(
    AUC = AUC,
    gdnss = gdnss,
    mu = mu,
    alpha = alpha
    ))
}

zetaIniOpt <- function(zeta, muIni, alphaIni, yi){
  yip <- Pz(muIni, alphaIni, zeta, Inf)
  return(abs(yi - yip))
}

CBMNLL <- function(mu, alpha, fi, ti){
  if (alpha > 1 || alpha < 0 ){
    return (-sum(.Machine$double.min.exp * c(fi, ti)))    
  }
  allParameters <- names(formals())
  zetaPos <- regexpr("zeta", allParameters)
  zeta <- unlist(mget(allParameters[which(zetaPos == 1)]))
  if (is.unsorted(zeta)){
    return (-sum(.Machine$double.min.exp * c(fi, ti)))
  }
  zeta <- c(-Inf, zeta, Inf)
  Q <- vector(length = length(fi))
  P <- Q
  for (i in 1:length(Q)){
    Q[i] <- Qz(zeta[i], zeta[i + 1])
    P[i] <- Pz(mu, alpha, zeta[i], zeta[i + 1])
  }

  L <- sum(log(c(Q, P)) * c(fi, ti))
  return(-L)
}

Qz <- function(zeta1, zeta2){
  return(pnorm(zeta2) - pnorm(zeta1))
}

Pz <- function(mu, alpha, zeta1, zeta2){
  P <- (1 - alpha) * (pnorm(zeta2) - pnorm(zeta1)) + alpha * (pnorm(zeta2, mean = mu) - pnorm(zeta1, mean = mu))
  return(P)
}

# pz <- function(mu, alpha, z){
#   p <- (1 - alpha) * dnorm(z) + alpha * dnorm(z, mean = mu)
#   return(p)
# }

CBMAucInt <- function(FPF, mu, alpha){
  z <- qnorm(FPF, lower.tail = FALSE)
  TPF <- Pz(mu, alpha, z, Inf)
  return(TPF)
}

GdnsFtCBM <- function(coef, fpBinnedTable, tpBinnedTable, K1, K2){
  mu <- coef[1]
  alpha <- coef[2]
  zeta <- c(-Inf, coef[3:length(coef)], Inf)
  FPF <- rep(NA, length(zeta) - 1)
  TPF <- FPF
  for (z in 2:length(zeta)){
    FPF[z - 1] <- Qz(zeta[z - 1], zeta[z])
    TPF[z - 1] <- Pz(mu, alpha, zeta[z - 1], zeta[z])
  }
  numFP <- FPF * K1
  numTP <- TPF * K2
  chiStat <- sum((numFP - fpBinnedTable)^2/numFP, (numTP - tpBinnedTable)^2/numTP)
  dgrfrdm <- length(fpBinnedTable) - 2
  pval <- 1 - pchisq(chiStat, dgrfrdm)
  return(pval)
}