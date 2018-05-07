Cov2Validation <- function(nBoots, z1ij, z2ij, DesiredNumBins, selI, aucEst = "RocFit")
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
  
  truth <- c(rep(0,K1),rep(1,K2)) 
  ret <- orCovariances(TRUE, Hratingsb, truth, aucEst)
  orgCov2 <- ret$cov2
  
  bsDir <- paste0("./REZ/",DataFile,"/","VAL/BS/")
  if (!dir.exists(bsDir)) dir.create(bsDir, recursive = TRUE)
  MyFile <- paste(paste0(bsDir ,"cov2BS"), "I", Istring, "J", Jstring, "BINS", DesiredNumBins, "BS", nBoots, aucEst, DataFile, sep = "_")
  if (!file.exists(MyFile)) { # as this takes a long time, let's save the results of the run
  registerDoMC(cores = detectCores())
  cov2 <- foreach (
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
      c(jb, cov$cov2, orgCov2)
    }
  save(file = MyFile, cov2)
  }else load(MyFile)
  return (cov2)
}
