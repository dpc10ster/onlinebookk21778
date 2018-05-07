CBMFitR <- function(dataset, selectedMod, selectedRdr){
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
  
  fp <- apply(NL, c(1, 2, 3), max)
  tp <- apply(LL, c(1, 2, 3), max)
  
  mu <- array(dim = c(I, J))
  alpha <- array(dim = c(I, J))
  AUC <- array(dim = c(I, J))
  gdnss <- array(dim = c(I, J))
  
  aucArray <- FigureOfMerit(dataset, "HrAuc")
  aucArray <- aucArray[selectedMod, selectedRdr]
  dim(aucArray) <- c(I, J)
  lenZetas <- length(unique(c(as.vector(fp), as.vector(tp)))) - 1
  zetas <- array(dim = c(I, J, lenZetas))
  for (i in 1:I){
    for (j in 1:J){
      scores <- sort(unique(c(fp[i, j, ], tp[i, j, ])))
      nBins <- length(scores)
      fpCounts <- rep(NA, nBins)
      tpCounts <- fpCounts
      for (b in 1:nBins){
        fpCounts[b] <- sum(fp == scores[b])
        tpCounts[b] <- sum(tp == scores[b])
      }
      
      FPF <- cumsum(rev(fpCounts)) / K1
      TPF <- cumsum(rev(tpCounts)) / K2
      
      alphaIni <- max(min(1 - (1 - rev(TPF)[which(rev(FPF) < 1)[1]]) / (1 - rev(FPF)[which(rev(FPF) < 1)[1]]), 0.99), 0.01)
      muIni <- qnorm(min((aucArray[i, j] - (1 - alphaIni)/2) / alphaIni, 0.99)) * sqrt(2)
      fpCountsCum <- cumsum(fpCounts)[1:(length(fpCounts) - 1)]
      fpCountsCum[fpCountsCum == 0] <- pnorm(minZeta)
      tpCountsCum <- cumsum(tpCounts)[1:(length(tpCounts) - 1)]
      tpCountsCum[tpCountsCum == 0] <- pnorm(minZeta, mean = muIni)
      
      zetaNorIni <- qnorm(fpCountsCum / K1)
      zetaAbnIni <- sapply(tpCountsCum / K2, IniZetaAbn, alpha = alphaIni, mu = muIni)
      for (z in 2:length(zetaNorIni)){
        zetaNorIni[z] <- min(max(zetaNorIni[z], zetaNorIni[z - 1] + 0.1), maxZeta)
        zetaAbnIni[z] <- min(max(zetaAbnIni[z], zetaAbnIni[z - 1] + 0.1), maxZeta)
      }
      for (z in (length(zetaNorIni) - 1):1){
        zetaNorIni[z] <- max(min(zetaNorIni[z], zetaNorIni[z + 1] - 0.1), minZeta)
        zetaAbnIni[z] <- max(min(zetaAbnIni[z], zetaAbnIni[z + 1] - 0.1), minZeta)
      }
      zetaIni <- (zetaNorIni + zetaAbnIni) / 2
      
      namesVector <- c("muFwd", "alphaFwd")
      for (t in 1:length(zetaIni))
        namesVector <- c(namesVector, paste0("zetaFwd", t))
      muIniFwd <- ForwardValue(muIni, minMu, maxMu)
      alphaIniFwd <- ForwardValue(alphaIni, minAlpha, maxAlpha)
      zetaIniFwd <- ForwardZetas(zetaIni)
      parameters <- c(list(muIniFwd, alphaIniFwd), as.list(zetaIniFwd))
      names(parameters) <- namesVector
      CBMNLLNew <- addArguments(CBMNLL, length(zetaIni))
      # ret <- mle2(CBMNLLTest, start = parameters, method = "BFGS", data = list(fi = fpCounts, ti = tpCounts))
      ret <- mle2(CBMNLLNew, start = parameters, method = "BFGS", data = list(fi = fpCounts, ti = tpCounts))
      mu[i, j] <- InverseValue(ret@coef[1], minMu, maxMu)
      alpha[i, j] <- InverseValue(ret@coef[2], minAlpha, maxAlpha)
      zetasTemp <- InverseZetas(ret@coef[3:length(ret@coef)])
      zetas[i, j, 1:length(zetasTemp)] <- zetasTemp
      # AUC[i, j] <- integrate(CBMAucInt, 0, 1, mu = mu[i, j], alpha = alpha[i, j])$value # XZ code 12/14/16
      # Browse[2]> integrate(CBMAucInt, 0, 1, mu = mu[i, j], alpha = alpha[i, j])$value
      # [1] 0.813925
      # Browse[2]> 0.5*(1-alpha[i, j]) + alpha[i, j]*pnorm(mu[i,j]/sqrt(2))
      # [1] 0.813925
      AUC[i, j] <- 0.5*(1-alpha[i, j]) + alpha[i, j]*pnorm(mu[i,j]/sqrt(2)) #dpc code 12/14/16
      gdnss[i, j] <- GdnsFtCBM(ret@coef, fpCounts, tpCounts, K1, K2)
    }
  }
  return(list(
    AUC = AUC,
    gdnss = gdnss,
    mu = mu,
    alpha = alpha,
    cutoffs = zetas
  ))
}

IniZetaAbn <- function(alpha, mu, tpCum){
  if (TPCumDiff(minZeta, alpha, mu, tpCum) >= 0){
    return(minZeta)
  }else if (TPCumDiff(maxZeta, alpha, mu, tpCum) <= 0){
    return(maxZeta)
  }else{
    return(uniroot(TPCumDiff, c(minZeta, maxZeta), alpha = alpha, mu = mu, tpCum = tpCum)$root)  
  }
}

TPCumDiff <- function(zeta, alpha, mu, tpCum){
  return((1 - alpha) * pnorm(zeta) + alpha * pnorm(zeta - mu) - tpCum)
}

CBMNLLTest <- function (muFwd, alphaFwd, fi, ti, zetaFwd1, zetaFwd2, zetaFwd3, 
                        zetaFwd4, zetaFwd5, zetaFwd6, zetaFwd7, zetaFwd8) 
{
  mu <- InverseValue(muFwd, minMu, maxMu)
  alpha <- InverseValue(alphaFwd, minAlpha, maxAlpha)
  allParameters <- names(formals())
  zetaPos <- regexpr("zeta", allParameters)
  zetaFwd <- unlist(mget(allParameters[which(zetaPos == 1)]))
  zeta <- InverseZetas(zetaFwd)
  zeta <- c(-Inf, zeta, Inf)
  Q <- vector(length = length(fi))
  P <- Q
  for (i in 1:length(Q)) {
    Q[i] <- Qz(zeta[i], zeta[i + 1])
    P[i] <- Pz(mu, alpha, zeta[i], zeta[i + 1])
  }
  L <- sum(log(c(Q, P)) * c(fi, ti))
  if (!is.finite(L)){
    dummyStop <- 1
  }
  return(-L)
}

CBMNLL <- function(muFwd, alphaFwd, fi, ti){
  mu <- InverseValue(muFwd, minMu, maxMu)
  alpha <- InverseValue(alphaFwd, minAlpha, maxAlpha)
  allParameters <- names(formals())
  zetaPos <- regexpr("zeta", allParameters)
  zetaFwd <- unlist(mget(allParameters[which(zetaPos == 1)]))
  zeta <- InverseZetas(zetaFwd)
  
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