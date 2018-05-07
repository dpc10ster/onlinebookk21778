LrocFoms <- function (zjk1, zjk2Cl, FPFValue) { 
  J <- length(zjk1[,1])
  PCL <- array(dim = J)
  ALroc <- array(dim = J)
  for (j in 1:J) {
    zk1 <- zjk1[j,]
    zk2Cl <- zjk2Cl[j,]    
    lroc <- LrocOperatingPointsFromRatings( zk1, zk2Cl )
    PCL[j] <- (approx(lroc$FPF, lroc$PCL, xout = FPFValue))$y
    ALroc[j] <- trapz(lroc$FPF, lroc$PCL)
  }
  return (list ( 
    PCL = PCL,
    ALroc = ALroc
  ))  
}