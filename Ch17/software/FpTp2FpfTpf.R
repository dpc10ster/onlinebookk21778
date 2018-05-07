FpTp2FpfTpf <- function (fp, tp) {
  zetas <- sort(unique(c(fp, tp)))
  nBins <- length(zetas)
  fpCounts <- rep(NA, nBins)
  tpCounts <- fpCounts
  for (bIndx in 1:nBins){
    fpCounts[bIndx] <- sum(fp == zetas[bIndx])
    tpCounts[bIndx] <- sum(tp == zetas[bIndx])
  }
  K1 <- length(fp)
  K2 <- length(tp)
  fpf <- cumsum(rev(fpCounts)) / K1
  tpf <- cumsum(rev(tpCounts)) / K2
  fpf <- fpf[-length(fpf)]
  tpf <- tpf[-length(tpf)]
  return(list(
    fpCounts = fpCounts,
    tpCounts = tpCounts,
    fpf <- fpf,
    tpf <- tpf,
    zetas = zetas
  ))
}

