simulateTheData <- function(simSelect, MyFile, S, I, J, K1, K2, ChiOrg, ChiOrgBar, 
                            mu, tau, ChiNew, varComp, DesiredNumBins, overWrite, isBinned = TRUE)
{
  truth <- c(rep(0,K1),rep(1,K2))
  Hratings <- array(dim = c(S,I,J,(K1+K2)))
  s <- 1
  if (overWrite || !file.exists(MyFile)) { # as this takes a long time, let's save the results of the run
    registerDoMC(cores = detectCores())
    saveSim <-  foreach (s = 1:S, .packages = c("lmf", "MASS"), .export = c("orCovariances")) %dorng% {
      fom <- array(dim = c(I, J))
      ret <- switch (simSelect,
                     CZ_MRMC = SimulateRocDataNew (I, J, K1, K2, ChiNew[s, , , , , ], DesiredNumBins),  
                     RM = RMSimulator(I, J, K1, K2, mu, tau, varComp, isBinned = TRUE, DesiredNumBins)
      )
      Hratings[s,,,1:K1] <- ret$z1b   
      Hratings[s,,,(K1+1):(K1+K2)] <- ret$z2b 
      HratingsTemp <- Hratings[s,,,]
      if (I == 1) {
        dim(HratingsTemp) <- c(I, dim(HratingsTemp))
      }
      saveSim <- orCovariances(TRUE, HratingsTemp, truth, aucEst)
      if (I == 1) {
        list(saveSim$cov2, saveSim$var, mean(saveSim$fomArray))
      } else {
        if (nhCondition == TRUE){
          c(list(saveSim$cov1, saveSim$cov2, saveSim$cov3, saveSim$var), as.list(apply(saveSim$fomArray, 1, mean)), list(saveSim$f, saveSim$ddf, saveSim$pValue))
        }else{
          c(list(saveSim$cov1, saveSim$cov2, saveSim$cov3, saveSim$var), as.list(apply(saveSim$fomArray, 1, mean)))
        }
      }
    }
    save(file = MyFile, saveSim)
  } else load(MyFile)
  
  return(saveSim)
}