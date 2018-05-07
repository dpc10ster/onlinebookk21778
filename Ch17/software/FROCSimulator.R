# changed input parameters to take lambda and BETA
FROCSimulator <- function(mu, lambda, nu, K1, K2, Nk22, zeta1){
  lambdaP <- lambda/mu
  nuP <- 1-exp(-nu*mu)
  nNL <- rpois(K1 + K2, lambdaP)
  maxNL <- max(nNL)
  NL <- array(dim = c(K1 + K2, maxNL))
  for (i in 1:(K1 + K2)){
    nl <- rnorm(nNL[i])
    nl <- nl[order(nl, decreasing = TRUE)]
    nl[nl < zeta1] <- NA
    NL[i, ] <- c(nl, rep(NA, maxNL - nNL[i]))
  }
  
  maxLes <- max(Nk22)
  LL <- array(dim = c(K2, maxLes))
  for (i in 1:K2){
    nLL <- rbinom(1, Nk22[i], nuP)
    ll <- rnorm(nLL, mu)
    ll <- ll[order(ll, decreasing = TRUE)]
    ll[ll < zeta1] <- NA
    LL[i, ] <- c(ll, rep(NA, maxLes - nLL))
  }
  
  lesID <- array(dim = c(K2, maxLes))
  lesWght <- array(dim = c(K2, maxLes))
  for (i in 1:K2){
    lesID[i, ] <- c(1:Nk22[i], rep(-Inf, maxLes - Nk22[i]))
    lesWght[i, ] <- c(rep(1 / Nk22[i], Nk22[i]), rep(-Inf, maxLes - Nk22[i]))
  }
  
  NL[is.na(NL)] <- -Inf
  LL[is.na(LL)] <- -Inf

  dim(NL) <- c(1, 1, K1 + K2, maxNL)
  dim(LL) <- c(1, 1, K2, maxLes)
  return(list(NL = NL, 
              LL = LL,
              lesionNum = Nk22,
              lesionID = lesID,
              lesionWeight = lesWght,
              dataType = "FROC",
              modalityID = "1",
              readerID = "1"
              ))
}