setupSimulator <- function (S, I, J, K1, K2, ChiOrgBar, SigmaChiOrgBar, orgData, CalSimuVals,
                            nhCondition, resampleChi, DesiredNumBins, RMVarStr, aucEst)
{
  varComp <-  NA
  ChiNew <- NA
  mu <- NA
  tau <- NA
  set.seed(seed)
  if (nhCondition == TRUE) nhStr <- "NH" else nhStr <- ""
  if (simSelect == "CZ_MRMC") {
    if (aucEst == "RocFit") 
    {
      retDir <- paste0("./REZ/", DataFile, "/VAL/RF/") 
    } else {
      retDir <- paste0("./REZ/", DataFile, "/VAL/TP/")
    }
  } else {
    retDir <- paste0("./REZ/", DataFile, "/RM/", RMVarStr, "/") 
  }
  retDir <- paste(paste0(retDir), "I", I, "J", J, "K1", K1, "K2", K2, sep = "_")
  retDir <- sub("_I", "I", retDir);retDir <- paste0(retDir,"/")
  if (!dir.exists(retDir)){
    dir.create(retDir, recursive = TRUE)
  }
  if (simSelect == "CZ_MRMC") {
      MyFile <- paste(paste0(retDir, nhStr), "BINS", DesiredNumBins, "S", S, "RSC", resampleChi, sep = "_")
      MyFile <- sub("_BINS", "BINS", MyFile)
      ChiNew <- SigmaChiNew (MyFile, S, I, J, ChiOrgBar, SigmaChiOrgBar) ## TEMP*DPC
  } else if (simSelect == "RM" ){
    MyFile <- paste(retDir, nhStr, "BINS", DesiredNumBins, "VarStr", RMVarStr, "S", S, "RSC", resampleChi, sep = "_")
    MyFile <- sub("_NH", "NH", MyFile);MyFile <- sub("__BINS", "BINS", MyFile)
    retRM <- SelectRMVarStr (I,CalSimuVals$AxOrgBar[5:6],RMVarStr,nhCondition)
    varComp  <- retRM$varComp
    mu <- retRM$mu
    tau <- retRM$tau
  } else {
    stop("Should not get here")
  }
  
  #   cat("Start ChiOrg dump\n") # this gets pasted to Excel sheet and unscrambled to get Tables 2 - 7
  #   done = array(0, dim = c(I,I,J,J))
  #   for (i1 in 1:I){
  #     for (j1 in 1:J){
  #       for (i2 in 1:I){
  #         for (j2 in 1:J){
  #           if (done[i1,i2,j1,j2] == 1) next
  #           cat(i1, j1, i2, j2, ChiNew[1,i1,i2,j1,j2,],"\n")
  #           done[i1,i2,j1,j2] <- 1
  #           done[i2,i1,j2,j1] <- 1
  #         }
  #       }
  #     }
  #   }
  #   cat("End ChiOrg dump\n\n")
  
  return(list(
    ChiNew = ChiNew,
    varComp = varComp,
    mu = mu,
    tau = tau,
    MyFile = MyFile
  ))
}