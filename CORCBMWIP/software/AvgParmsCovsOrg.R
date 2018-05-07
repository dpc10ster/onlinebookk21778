AvgParmsCovsOrg <- function(parms,Cov,Ax,Ay,RhoAxAy,done)
{
  I <- length(parms[,1,1,1,1])   
  J <- length(parms[1,1,,1,1])
  dim1 <- 6 # renumbering; keeping modalities separate; current code is limited to I = 2 (TEMP*DPC)
  ChiBar <- array(0, dim = c(dim1,6))
  SigmaChiBar <- array(0,dim=c(dim1,6,6))
  AxAvg <- array(0, dim = dim1)
  AyAvg <- array(0, dim = dim1)
  RhoAxAyAvg <- array(0, dim = dim1)
  AvgCount <- array(0, dim = dim1)
  
  RhoAxAyAvg[5:6] <- 1
  for (i1 in 1:I) {
    for (i2 in 1:I) {        
      for (j1 in 1:J) {
        for (j2 in 1:J) {
          if (done[i1,i2,j1,j2] == 0) next
          p <- PairingType(i1, i2, j1, j2)
          if (p < 5) {
            ChiBar[p,] <- ChiBar[p,] + parms[i1,i2,j1,j2,]
            SigmaChiBar[p,,] <- SigmaChiBar[p,,] + Cov[i1,i2,j1,j2,,]
            AxAvg[p] <- AxAvg[p] + Ax[i1,i2,j1,j2]
            AyAvg[p] <- AyAvg[p] + Ay[i1,i2,j1,j2]
            RhoAxAyAvg[p] <- RhoAxAyAvg[p] + RhoAxAy[i1,i2,j1,j2]
            AvgCount[p] <- AvgCount[p] + 1
          } else {
            ChiBar[p,] <- ChiBar[p,] + parms[i1,i2,j1,j2,]
            SigmaChiBar[p,1:2,1:2] <- SigmaChiBar[p,1:2,1:2] + Cov[i1,i2,j1,j2,1:2,1:2]
            AxAvg[p] <- AxAvg[p] + Ax[i1,i2,j1,j2]
            AyAvg[p] <- AxAvg[p]
            AvgCount[p] <- AvgCount[p] + 1
          }
        }
      }
    }
  }
  
  for (p in 1:(4+I)) {
    ChiBar[p,] <- ChiBar[p,]/AvgCount[p]
    SigmaChiBar[p,,] <- SigmaChiBar[p,,]/AvgCount[p]
    AxAvg[p] <- AxAvg[p]/AvgCount[p]
    AyAvg[p] <- AyAvg[p]/AvgCount[p]
    if (p < 5) RhoAxAyAvg[p] <- RhoAxAyAvg[p]/AvgCount[p]
  }

  SigmaChiBar[5:6,c(3,4),c(3,4)] <- SigmaChiBar[5:6,c(1,2),c(1,2)]
  SigmaChiBar[5:6,5,5] <- 1
  SigmaChiBar[5:6,6,6] <- 1
  
  return(list(
    ChiBar = ChiBar, 
    SigmaChiBar= SigmaChiBar,
    AxAvg = AxAvg,
    AyAvg = AyAvg,
    RhoAxAyAvg = RhoAxAyAvg
  ))
  
}