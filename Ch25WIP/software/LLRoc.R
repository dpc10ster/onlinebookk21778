oLoopVaryLambdaRoc <- function(lambdaPFwd, fb, tb, K1, K2, nLesDistr, AUCCbm, nuPIni, zetaIni, lambdaPLower, iniCal){  
  
  lambdaP <- InverseLamdba(lambdaPFwd, lambdaPLower)
  
  ret <- iLoopVaryNuZetasRoc(lambdaP, fb, tb, K1, K2, nLesDistr, AUCCbm, nuPIni, zetaIni, iniCal)
  
  return(ret$mleRet@min)
}

iLoopVaryNuZetasRoc <- function(lambdaP, fb, tb, K1, K2, nLesDistr, AUCCbm, nuPIni, zetaIni, iniCal){
  minNuP <- RJafrocEnv$minNuP
  maxNuP <- RJafrocEnv$maxNuP
  
  maxMu <- RJafrocEnv$maxMu
  minMu <- RJafrocEnv$minMu
  
  # Look for better lower bound nuPLower using maxMu. 
  nuPLower <- uniroot(DiffAucSmMinusAucCbm, c(minNuP, maxNuP), mu = maxMu, lambdaP = lambdaP, 
                      nLesDistr = nLesDistr, 
                      AUCCbm = AUCCbm, 
                      tol = .Machine$double.eps^0.5)$root
  
  # Above implies that if nuP(true) is lower than nuPLower, it is impossible to make AUCSm equal to AUCCbm
  
  if (DiffAucSmMinusAucCbm(minMu, lambdaP, maxNuP, nLesDistr, AUCCbm) > 0){ 
    # Above implies AUCSm > AUCCbm even using minMu.
    # This means maxNuP is too large. 
    # Need to find a better (i.e., smaller) upper bound nuPUpper, again using minMu
    nuPUpper <- uniroot(DiffAucSmMinusAucCbm, c(nuPLower, maxNuP), mu = minMu, lambdaP = lambdaP, 
                        nLesDistr = nLesDistr, AUCCbm = AUCCbm, tol = .Machine$double.eps^0.5)$root
  }else{
    # assign default upper bound
    nuPUpper <- maxNuP
  }
  
  if (nuPUpper - nuPIni <= .Machine$double.eps^0.5){
    # nuPIni is too close (within tol), or larger, than nuPUpper, 
    # then predicted AUC(nuPIni) is slightly smaller or 
    # larger than predicted AUC(nuPUpper) 
    # but still larger than AUCCbm. Then mu cannot be found.
    # so we reduced nuPIni to a safe value below nuPUpper
    nuPIni <- nuPUpper - 2 * .Machine$double.eps^0.5
  }else if (nuPIni - nuPLower <= .Machine$double.eps^0.5){
    # Similar argument, but applies to lower bound.
    nuPIni <- nuPLower + 2 * .Machine$double.eps^0.5
  }
  
  # zetaIni <- rep(NA, length(fb) - 1) 
  # revFPF <- rev(cumsum(rev(fb)) / K1)[-1] # gets rid of the 1, leaving the highest non-trivial FPF value
  # minZeta <- RJafrocEnv$minZeta
  # maxZeta <- RJafrocEnv$maxZeta
  # 
  # # following follows from expression for FPF and relation between erf and Phi functions
  # phiZeta <- (log(1 - revFPF)/lambdaP + 1) # xz was right
  # for (z in 2:length(phiZeta)){
  #   phiZeta[z] <- max(phiZeta[z - 1] + 0.01, phiZeta[z])
  # } 
  # if (any(phiZeta < 0.01)) {
  #   # dpc!! if this error occurs, the trial lambda is too small
  #   for (z in 1:length(phiZeta)){
  #     if (phiZeta[z] < 0.01){
  #       # impossible to found Gaussian quantile
  #       if (z == 1){
  #         # assign a small probability
  #         phiZeta[z] <- 0.01
  #       }else{
  #         # assign a small probability and larger than the previous cutoff
  #         phiZeta[z] <- phiZeta[z - 1] + 0.01
  #       }
  #     }else{
  #       # make sure every cutoff is larger than the previous one
  #       phiZeta[z] <- max(phiZeta[z - 1] + 0.01, phiZeta[z])
  #     }
  #   }
  # }
  # 
  # if (any(phiZeta > 0.99)){
  #   for (z in length(phiZeta):1){
  #     if (phiZeta[z] > 0.99){
  #       if (z == length(phiZeta)){
  #         phiZeta[z] <- 0.99
  #       }else{
  #         phiZeta[z] <- phiZeta[z + 1] - 0.01
  #       }
  #     }else{
  #       phiZeta[z] <- min(phiZeta[z + 1] - 0.01, phiZeta[z])
  #     }
  #   }
  # }
  # zetaIni <- qnorm(phiZeta) 
  zetaIniFwd <- ForwardZetas(zetaIni)
  
  nuPIniFwd <- ForwardNu(nuPIni, nuPLower, nuPUpper)
  namesVector <- NULL
  for (t in 1:length(zetaIniFwd))
    namesVector <- c(namesVector, paste0("zetaFwd", t))
  namesVector <- c("nuPFwd", namesVector)
  parameters <- c(nuPIniFwd, zetaIniFwd)
  names(parameters) <- namesVector
  
  LLRocNew <- addArguments(LLRoc, length(zetaIniFwd))
  names(parameters) <- namesVector
  if (iniCal){
    logLikIni <- LLRocNew(lambdaP, nuPIniFwd, fb, tb, nLesDistr, AUCCbm, nuPLower, nuPUpper, zetaIni)
    return(logLikIni = logLikIni)
  }else{
    # now vary the zetas 
    if (any(!is.finite(parameters))){
      dummyStop <- 1
    }
    # ret <- mle2(minuslogl = LLRocTest, start = as.list(parameters), method = "BFGS", 
    #             data = list(lambdaP = lambdaP, fb = fb, tb = tb, nLesDistr = nLesDistr,
    #                         AUCCbm = AUCCbm, nuPLower = nuPLower, nuPUpper = nuPUpper) )
    ret <- mle2(minuslogl = LLRocNew, start = as.list(parameters), method = "BFGS",
                data = list(lambdaP = lambdaP, fb = fb, tb = tb, nLesDistr = nLesDistr,
                            AUCCbm = AUCCbm, nuPLower = nuPLower, nuPUpper = nuPUpper) )
    cat("logLikInner = ", ret@min, ", lambdaP = ", lambdaP, "coef = ", ret@fullcoef, "\n\n\n")
    if (DEBUG) cat("logLikInner = ", ret@min, ", lambdaP = ", lambdaP, "\n")
    return(list(
      mleRet = ret,
      nuPLower = nuPLower,
      nuPUpper = nuPUpper
    ))
  }
}

# this function calculates -LL for specified zetas
LLRoc <- function(lambdaP, nuPFwd, fb, tb, nLesDistr, AUCCbm, nuPLower, nuPUpper, zetaIni){
  maxMu <- RJafrocEnv$maxMu
  minMu <- RJafrocEnv$minMu
  if (missing(zetaIni)){
    allParameters <- names(formals())
    zetaPos <- grep("zetaFwd", allParameters)
    zeta <- unlist(mget(allParameters[zetaPos]))
    zeta <- InverseZetas(zeta)
  }else{
    zeta <- zetaIni
  }
  nuP <- InverseNu(nuPFwd, nuPLower, nuPUpper)
  # if (DiffAucSmMinusAucCbm(minMu, lambdaP, nuP, nLesDistr, AUCCbm) * DiffAucSmMinusAucCbm(maxMu, lambdaP, nuP, nLesDistr, AUCCbm) > 0){
  #   dummyStop <- 1
  # }
  
  # Same problem for nuPIni
  if (nuPUpper - nuP <= .Machine$double.eps^0.5){
    nuP <- nuPUpper - 2 * .Machine$double.eps^0.5
  }else if (nuP - nuPLower <= .Machine$double.eps^0.5){
    nuP <- nuPLower + 2 * .Machine$double.eps^0.5
  }
  mu <- uniroot(DiffAucSmMinusAucCbm, c(minMu, maxMu), lambdaP = lambdaP, nuP = nuP,
                nLesDistr = nLesDistr, AUCCbm = AUCCbm, tol = .Machine$double.eps^0.5)$root
  FPF <- xROC(zeta, lambdaP)
  TPF <- yROC(zeta, mu, lambdaP, nuP, nLesDistr)
  FPFTerms <- c(1, FPF) - c(FPF, 0)
  TPFTerms <- c(1, TPF) - c(TPF, 0)
  FPFTerms[ FPFTerms < 1e-15 ] <-  1e-15
  TPFTerms[ TPFTerms < 1e-15 ] <-  1e-15
  L <- sum(log(c(FPFTerms, TPFTerms)) * c((fb), (tb))) 
  cat(L, nuP, zeta, mu, "\n")
  return (-L)
}

LLRocTest <- function (lambdaP, nuPFwd, fb, tb, nLesDistr, AUCCbm, nuPLower, 
          nuPUpper, zetaIni, zetaFwd1, zetaFwd2, zetaFwd3, zetaFwd4, 
          zetaFwd5) 
{
  minMu <- RJafrocEnv$minMu
  if (missing(zetaIni)) {
    allParameters <- names(formals())
    zetaPos <- grep("zetaFwd", allParameters)
    zeta <- unlist(mget(allParameters[zetaPos]))
    zeta <- InverseZetas(zeta)
  }
  else {
    zeta <- zetaIni
  }
  nuP <- InverseNu(nuPFwd, nuPLower, nuPUpper)
  if (nuPUpper - nuP <= .Machine$double.eps^0.5) {
    nuP <- nuPUpper - 2 * .Machine$double.eps^0.5
  }
  else if (nuP - nuPLower <= .Machine$double.eps^0.5) {
    nuP <- nuPLower + 2 * .Machine$double.eps^0.5
  }
  mu <- uniroot(DiffAucSmMinusAucCbm, c(minMu, maxMu), lambdaP = lambdaP, 
                nuP = nuP, nLesDistr = nLesDistr, AUCCbm = AUCCbm, tol = .Machine$double.eps^0.5)$root
  FPF <- xROC(zeta, lambdaP)
  TPF <- yROC(zeta, mu, lambdaP, nuP, nLesDistr)
  FPFTerms <- c(1, FPF) - c(FPF, 0)
  TPFTerms <- c(1, TPF) - c(TPF, 0)
  FPFTerms[FPFTerms < 1e-15] <- 1e-15
  TPFTerms[TPFTerms < 1e-15] <- 1e-15
  L <- sum(log(c(FPFTerms, TPFTerms)) * c((fb), (tb)))
  if (!is.finite(L)){
    dummyStop <- 1
  }
  return(-L)
}