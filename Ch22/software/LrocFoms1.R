LrocFoms <- function (zk1, zk2Cl, FPFValue) { 
  lroc <- LrocOperatingPointsFromRatings( zk1, zk2Cl )
  PCL <- (approx(lroc$FPF, lroc$PCL, xout = FPFValue))$y # computes PCL @ FPFValue
  tempFpf <-c(lroc$FPF[lroc$FPF < FPFValue],FPFValue)
  tempPcl <-c(lroc$PCL[lroc$FPF < FPFValue],PCL)
  ALroc <- trapz(tempFpf, tempPcl) # computes trapezoidal area under LROC (0 to FPFValue)
  return (list ( 
    PCL = PCL,
    ALroc = ALroc
  ))  
}