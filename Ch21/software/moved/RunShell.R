RunShell <- function(DesiredNumBins,J,K1,K2,z1,z2)
{
  
  # matrix order
  # COMMON/BLK7/AX,BX,AY,BY,R,RHO,T(11),U(11)
  DIM <- 6+2*(DesiredNumBins-1)#ax,ay,bx,by,rhon,rhos,ctffsx[], ctffsy[] 
  CovArr <- array(dim=c(J,DIM,DIM)) 
  ax <- array(dim=J);ay <- array(dim=J);bx <- array(dim=J);by <- array(dim=J)
  Ax <- array(dim=J);Ay <- array(dim=J);StdAx <- array(dim=J);StdAy <- array(dim=J);RhoAxAy <- array(dim=J)
  rhon <- array(dim=J);rhos <- array(dim=J)
  #zetan <- array(dim=c(J,4));zetas <- array(dim=c(J,4))
  done <- array(dim=J) 
  for (j in 1:J) {
    OutLines  <- array("", dim = K1+6)
    OutLines[1] <- "TYPICAL EXAMPLE OF 5-CATEGORY DATA --X"
    OutLines[2] <- max(z1,z2) 
    OutLines[3] <- max(z1,z2) 
    OutLines[4] <- "TYPICAL EXAMPLE OF 5-CATEGORY DATA --Y"
    OutLines[5] <- max(z1,z2) 
    OutLines[6] <- max(z1,z2) 
    for (k1 in 1:K1) {
      OutLines[6+k1] <- paste(z1[1,j,k1], z1[2,j,k1], sep =" ")
    }
    OutLines[6+k1+1] <- "*"
    for (k2 in 1:K2) {
      OutLines[6+k1+1+k2] <- paste(z2[1,j,k2], z2[2,j,k2], sep =" ")
    }
    OutLines[6+k1+1+k2+1] <- "*"
    OutLines[6+k1+1+k2+2] <- "area"
    
    writeLines(OutLines,"./CorROCII/Debug/1R2MData.txt")
    if (.Platform$OS.type != "windows") {
      system("./CorROCII.sh > NULL")
    } else {
      system("./CorROCII/Debug/CORROC2.bat > NULL")
    }
    unlink(NULL)
    filename <- "./CorROCII/Debug/CorrocIIoutput.txt"    
    ret <- ReadCorROCII(filename)
    
    if (length(ret) == 1) done[j] <- 0 else {
      done[j] <- 1     
      ax[j] <- ret$ax
      bx[j] <- ret$bx   
      ay[j] <- ret$ay
      by[j] <- ret$by
      rhon[j] <- ret$rhon
      rhos[j] <- ret$rhos
      #zetan[j,1:length(ret$zetan)] <- ret$zetan
      #zetas[j,1:length(ret$zetas)] <- ret$zetas
      Cov <- ret$Cov
      CovArr[j,1:length(Cov[1,]),1:length(Cov[1,])] <- Cov 
      Ax[j] <- ret$Ax
      Ay[j] <- ret$Ay
      StdAx[j] <- ret$StdAx
      StdAy[j] <- ret$StdAy
      RhoAxAy[j] <- ret$RhoAxAy 
      pValue  <- ret$pValue
      #cat("readers = ",j,"\n")
    }
  }
  
  if (all(done == 0)) return(-1)
  ax <- ax[which(done == 1)]
  bx <- bx[which(done == 1)]   
  ay <- ay[which(done == 1)]
  by <- by[which(done == 1)]
  rhon <- rhon[which(done == 1)]
  rhos <- rhos[which(done == 1)] 
  Ax <- Ax[which(done == 1)]
  Ay <- Ay[which(done == 1)]
  StdAx <- StdAx[which(done == 1)]
  StdAy <- StdAy[which(done == 1)]
  RhoAxAy <- RhoAxAy[which(done == 1)]
  CovArr <- CovArr[which(done == 1),,]
  
  TrtMeansOrig <- c(mean(Ax),mean(Ay))
  
  if (J > 1) CovOrig  <- apply(CovArr[,1:6,1:6], c(2,3), mean) else CovOrig <- Cov[1:6,1:6]
  axOrig <- mean(ax)
  bxOrig <- mean(bx)
  ayOrig <- mean(ay)
  byOrig <- mean(by)
  rhonOrig <- mean(rhon)
  rhosOrig <- mean(rhos)
  AxOrig <- mean(Ax);AyOrig <- mean(Ay)
  StdAxOrig <- mean(StdAx);StdAyOrig <- mean(StdAy)
  RhoAxAyOrig <- mean(RhoAxAy)
  
  #   cat("Original data values:, \n")
  #   cat("axOrig=", axOrig, "\n") 
  #   cat("bxOrig=", bxOrig, "\n") 
  #   cat("ayOrig=", ayOrig, "\n") 
  #   cat("byOrig=", byOrig, "\n") 
  #   cat("rhonOrig=", rhonOrig, "\n") 
  #   cat("rhosOrig=", rhosOrig, "\n") 
  #   cat("AreaOrig=", mean(c(AxOrig, AyOrig)), "\n") 
  
  return(list(done = done, 
              axOrig = axOrig,
              bxOrig = bxOrig,
              ayOrig = ayOrig,
              byOrig = byOrig,
              rhonOrig = rhonOrig,
              rhosOrig = rhosOrig,
              CovOrig = CovOrig,
              pValue = pValue))
  
}