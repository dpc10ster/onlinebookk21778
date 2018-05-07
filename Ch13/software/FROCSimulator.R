# changed input parameters to take lambda and BETA
FROCSimulator <- function(mu, lambda, nu, K1, K2, Lk2, zeta1){
  lambdaP <- lambda/mu
  nuP <- 1-exp(-nu*mu)
  nNL <- rpois(K1 + K2, lambdaP)
  maxNL <- max(nNL)
  NL <- array(dim = c(K1 + K2, maxNL))
  for (k in 1:(K1 + K2)){
    znl <- rnorm(nNL[k])
    znl <- znl[order(znl, decreasing = TRUE)]
    znl[znl < zeta1] <- NA
    NL[k, ] <- c(znl, rep(NA, maxNL - nNL[k]))
  }
  
  maxLes <- max(Lk2)
  LL <- array(dim = c(K2, maxLes))
  for (k2 in 1:K2){
    nLL <- rbinom(1, Lk2[k2], nuP)
    zll <- rnorm(nLL, mu)
    zll <- zll[order(zll, decreasing = TRUE)]
    zll[zll < zeta1] <- NA
    LL[k2, ] <- c(zll, rep(NA, maxLes - nLL))
  }
  
  lesID <- array(dim = c(K2, maxLes))
  lesWght <- array(dim = c(K2, maxLes))
  for (k2 in 1:K2){
    lesID[k2, ] <- c(1:Lk2[k2], rep(-Inf, maxLes - Lk2[k2]))
    lesWght[k2, ] <- c(rep(1 / Lk2[k2], Lk2[k2]), rep(-Inf, maxLes - Lk2[k2]))
  }
  
  NL[is.na(NL)] <- -Inf
  LL[is.na(LL)] <- -Inf

  dim(NL) <- c(1, 1, K1 + K2, maxNL)
  dim(LL) <- c(1, 1, K2, maxLes)
  return(list(NL = NL, 
              LL = LL,
              lesionNum = Lk2,
              lesionID = lesID,
              lesionWeight = lesWght,
              dataType = "FROC",
              modalityID = "1",
              readerID = "1"
              ))
}