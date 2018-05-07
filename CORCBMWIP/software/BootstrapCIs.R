BootstrapCIs <- function(nBoots, z1ij, z2ij, DesiredNumBins, selI, aucEst = "RocFit")
{  
  I <- dim(z1ij)[1];J <- dim(z1ij)[2];K1 <- dim(z1ij)[3];K2 <- dim(z2ij)[3]
  Hratingsb <- array(dim=c(I,J,(K1+K2)))
  for (i in 1:I){
    for (j in 1:J){
      ret <- ToIntegerRatings(z1ij[i, j, ],z2ij[i, j, ],DesiredNumBins)
      Hratingsb[i,j,1:K1] <- z1ij[i, j, ]   
      Hratingsb[i,j,(K1+1):(K1+K2)] <- z2ij[i, j, ]
      Hratingsb[i,j,1:K1] <- ret$f[1:K1]   
      Hratingsb[i,j,(K1+1):(K1+K2)] <- ret$t[1:K2]
    }
  }
  
  # no need to use overwrite flag here
  truth <- c(rep(0,K1),rep(1,K2)) 
  bsDir <- paste0("./REZ/",DataFile,"/","VAL/BS/")
  if (!dir.exists(bsDir)) dir.create(bsDir, recursive = TRUE)
  MyFile <- paste(paste0(bsDir ,"varCovBS"), "I", Istring, "J", Jstring, "BINS", DesiredNumBins, "BS", nBoots, aucEst, DataFile, sep = "_")
  if (!file.exists(MyFile)) { # as this takes a long time, let's save the results of the run
    registerDoMC(cores = detectCores())
    varComp <- foreach (
      b = 1:nBoots, 
      .combine = "rbind", 
      .export = c("orCovariances")) %dorng%  {
        kb1 <- ceiling(runif(K1) * K1) 
        kb2 <- ceiling(runif(K2) * K2) + K1
        jb <- ceiling(runif(J) * J)
        HratingsbB <- Hratingsb[ , jb, c(kb1, kb2)]
        if (length(selI) == 1)
          dim(HratingsbB) <- c(1, dim(HratingsbB))
        cov <- orCovariances(TRUE, HratingsbB, truth, aucEst)
        if (length(selI) == 1){
          as.numeric(c(cov$cov2, cov$var, mean(cov$fomArray)))
        }else{
          as.numeric(c(cov$cov1, cov$cov2, cov$cov3, cov$var, apply(cov$fomArray, 1, mean)))
        }
      }
    CI <- array(dim = c(ncol(varComp), 2))
    bsMean <- rep(NA, ncol(varComp))
    bsVar <- rep(NA, ncol(varComp))
    for (n in 1:ncol(varComp)) {
      CI[n, ] <- quantile(varComp[, n], c(0.025, 0.975))
      bsMean[n] <- mean(varComp[, n])
      bsVar[n] <- var(varComp[, n])
    }
    
    ret <- orCovariances(TRUE, Hratingsb, truth, aucEst)
    if (length(selI) == 1){
      ret <- as.numeric(c(ret$cov2, ret$var, mean(ret$fomArray[1,]), mean(ret$fomArray)))
    }else{
      ret <- as.numeric(c(ret$cov1, ret$cov2, ret$cov3, ret$var, apply(ret$fomArray, 1, mean)))
    }
    save(file = MyFile, varComp, CI, bsMean, bsVar, ret)
  } else load(MyFile)
  
  return (list(
    I=length(selI),
    J=J,
    K1=K1,
    K2=K2,
    orgVarCompFom = ret,
    CI = CI, 
    bsMean = bsMean,
    bsVar = bsVar))
}
