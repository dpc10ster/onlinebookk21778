CreateNHParameters <- function(calI, calJ, ChiOrgBar, SigmaChiOrgBar)
{
  
  ChiOrgBarNH <- ChiOrgBar * 0
  SigmaChiOrgBarNH <- SigmaChiOrgBar * 0
  
  for (p in 1:2)
  {  
    ChiOrgBarNH[p,c(1,3)] <- (ChiOrgBar[p,1]+ChiOrgBar[p,3])/2
    ChiOrgBarNH[p,c(2,4)] <- (ChiOrgBar[p,2]+ChiOrgBar[p,4])/2
    ChiOrgBarNH[p,c(5,6)] <- ChiOrgBar[p,c(5,6)]
  }
  p <- 3
  ChiOrgBarNH[p,c(1,3)] <- (ChiOrgBar[p,1] + ChiOrgBar[p,3]+ChiOrgBar[p+1,1] + ChiOrgBar[p+1,3])/4
  ChiOrgBarNH[p,c(2,4)] <- (ChiOrgBar[p,2] + ChiOrgBar[p,4]+ChiOrgBar[p+1,2] + ChiOrgBar[p+1,4])/4
  ChiOrgBarNH[p,5] <- (ChiOrgBar[p,5]+ChiOrgBar[p+1,5])/2
  ChiOrgBarNH[p,6] <- (ChiOrgBar[p,6]+ChiOrgBar[p+1,6])/2
  ChiOrgBarNH[4,] <- ChiOrgBarNH[3,]
  p <- 5
  ChiOrgBarNH[p,c(1,3)] <- (ChiOrgBar[p,1] + ChiOrgBar[p,3]+ChiOrgBar[p+1,1] + ChiOrgBar[p+1,3])/4
  ChiOrgBarNH[p,c(2,4)] <- (ChiOrgBar[p,2] + ChiOrgBar[p,4]+ChiOrgBar[p+1,2] + ChiOrgBar[p+1,4])/4
  ChiOrgBarNH[p,c(5,6)] <- ChiOrgBar[p,c(5,6)]
  ChiOrgBarNH[6,] <- ChiOrgBarNH[5,]
 
  SigmaChiOrgBarNH[1:2,,] <- SigmaChiOrgBar[1:2,,] 
  SigmaChiOrgBarNH[3,,] <- (SigmaChiOrgBar[3,,]+SigmaChiOrgBar[4,,])/2
  SigmaChiOrgBarNH[4,,] <- SigmaChiOrgBarNH[3,,]
  SigmaChiOrgBarNH[5,,] <- (SigmaChiOrgBar[5,,]+SigmaChiOrgBar[6,,])/2
  SigmaChiOrgBarNH[6,,] <- SigmaChiOrgBarNH[5,,]
  
  ChiOrgBar <- ChiOrgBarNH
  SigmaChiOrgBar <- SigmaChiOrgBarNH
  
  return(list(
    ChiOrgBar = ChiOrgBar,
    SigmaChiOrgBar = SigmaChiOrgBar
  ))
}