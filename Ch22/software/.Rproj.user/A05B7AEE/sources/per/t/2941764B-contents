# NOTE: two versions yield identical results
# older version
LrocOperatingPointsFromRatings1 <- function( zk1, zk2Cl ) 
{
  FPF <- 1
  PCL <- NULL
  zk11 <- zk1;zk2Cl1 <- zk2Cl
  while(1) {
    cutoff <- min( c( zk11, zk2Cl1 ) )
    zk11 <- zk1[ zk1 > cutoff ]
    zk2Cl1 <- zk2Cl[ zk2Cl > cutoff ]
    FPF1 <- length( zk11 ) / length( zk1 )
    PCL1 <- length( zk2Cl1 ) / length( zk2Cl )
    FPF <- c( FPF, FPF1 )
    if (length(PCL) == 0) PCL <- c(PCL1, PCL1) else PCL <- c( PCL, PCL1 )
    if( FPF1 == 0 && PCL1 == 0 ) {
      break
    }
  }  
  return( list(
    FPF = FPF[length(FPF):1],
    PCL = PCL[length(PCL):1]
  ) )  
}


# cleaner version
LrocOperatingPointsFromRatings <- function( zk1, zk2Cl ) 
{
  bins <- sort(unique(c(zk1,zk2Cl)))
  nBins <- length(bins)
  
  fpCounts <- array(0, dim = nBins)
  clCounts <- array(0, dim = nBins)
  
  for (b in 1:nBins){
    fpCounts[b] <- sum(zk1 == bins[b])
    clCounts[b] <- sum(zk2Cl == bins[b])
  }
  
  FPF <- cumsum(rev(fpCounts)) / length(zk1)
  PCL <- cumsum(rev(clCounts)) / length(zk2Cl)
  FPF <- FPF[-length(FPF)]
  PCL <- PCL[-length(PCL)]
  FPF <- c(0, FPF, 1) # add origin and largest value
  PCL <- c(0, PCL, PCL[length(PCL)]) # add origin and max PCL value
  return( list(
    FPF = FPF,
    PCL = PCL
  ) )  
}