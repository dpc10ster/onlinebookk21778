source("ReadCorROCII.R")
RunCorrocOnPairedData <- function(z1b,z2b,DesiredNumBins) 
{
  K1 <- length(z1b[1,]);K2 <- length(z2b[1,])
  OutLines  <- array("", dim = K1+6)
  OutLines[1] <- "TYPICAL EXAMPLE OF 5-CATEGORY DATA --X"
  OutLines[2] <- max(z1b,z2b) 
  OutLines[3] <- max(z1b,z2b) 
  OutLines[4] <- "TYPICAL EXAMPLE OF 5-CATEGORY DATA --Y"
  OutLines[5] <- max(z1b,z2b) 
  OutLines[6] <- max(z1b,z2b) 
  for (k1 in 1:K1) {
    OutLines[6+k1] <- paste(z1b[1,k1], z1b[2,k1], sep =" ")
  }
  OutLines[6+k1+1] <- "*"
  for (k2 in 1:K2) {
    OutLines[6+k1+1+k2] <- paste(z2b[1,k2], z2b[2,k2], sep =" ")
  }
  OutLines[6+k1+1+k2+1] <- "*"
  OutLines[6+k1+1+k2+2] <- "area"
  
  filenameOut <- "./CorROCII/Debug/CorrocIIoutput.txt"
  filenameInp <- "./CorROCII/Debug/1R2MData.txt"  
  file.remove(filenameInp)     
  writeLines(OutLines,filenameInp)   
  if (Sys.info()["sysname"] == "Windows"){
    system2("CorROCII/Debug/CORROC2.bat")
  }else{
    system2("./CorROCII.sh", stdout = NULL) 
    # without the second argument the program stops here; 12/31/17
  }
  ret <- ReadCorROCII(filenameOut)
  parms <- array(dim = 6)
  Cov <- array(dim=c(6,6))
  if (length(ret) == 1) {
    done <- 0
    Ax <-  NA
    Ay <-  NA
    RhoAxAy <-  NA
  } else {
    done <- 1 
    parms[1] <- ret$ax
    parms[2] <- ret$bx   
    parms[3] <- ret$ay
    parms[4] <- ret$by
    parms[5] <- ret$rhon
    parms[6] <- ret$rhos
    Cov <- ret$Cov[1:6,1:6]
    Ax <- ret$Ax
    Ay <- ret$Ay
    pValue  <- ret$pValue
    RhoAxAy <- ret$RhoAxAy
  }
  cat("\n\n")
  return(list( 
    done = done,
    parms = parms,
    Cov = Cov,
    Ax = Ax,
    Ay = Ay,
    pValue = pValue,
    RhoAxAy=RhoAxAy
  ))
  
}

