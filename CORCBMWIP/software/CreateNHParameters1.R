# CreateNHParameters.R
CreateNHParameters <- function(calI, calJ, ChiOrgBar, SigmaChiOrgBar)
{
  ChiOrgBarNH = array(dim = dim(ChiOrgBar))
  ChiOrgBar <- ChiOrgBarNH
  for (p in (calI+3):(calI+4))
  {  
    ChiOrgBarNH[p,c(1,2,3,4)] <- (ChiOrgBar[p,c(1,2)] + ChiOrgBar[p,c(3,4)])/2
    ChiOrgBarNH[p,c(5,6)] <- ChiOrgBar[p,c(5,6)]
    ChiOrgBar[p,] <- (ChiOrgBarNH[(calI+3),]+ChiOrgBarNH[(calI+4),])/2
  }
  
  SigmaChiOrgBarNH <- SigmaChiOrgBar
  for (i1 in 1:calI) {
    for (i2 in 1:calI) {        
      for (j1 in 1:calJ) {
        for (j2 in 1:calJ) {
          p <- PairingType(i1, i2, j1, j2)
          SigmaChiOrgBarNH[p,1,1] <- (SigmaChiOrgBar[p,1,1] + SigmaChiOrgBar[p,3,3])/2
          SigmaChiOrgBarNH[p,3,3] <- SigmaChiOrgBarNH[p,1,1]
          SigmaChiOrgBarNH[p,2,2] <- (SigmaChiOrgBar[p,2,2] + SigmaChiOrgBar[p,4,4])/2
          SigmaChiOrgBarNH[p,4,4] <- SigmaChiOrgBarNH[p,2,2]
        }
      }
    }
  }
  SigmaChiOrgBar <- SigmaChiOrgBarNH
  for (p in (calI+1):(calI+2))
  {  
    SigmaChiOrgBar[p,,] <- (SigmaChiOrgBarNH[(calI+1),,]+SigmaChiOrgBarNH[(calI+2),,])/2
  }
  
  return(list(
    ChiOrgBar = ChiOrgBar,
    SigmaChiOrgBar = SigmaChiOrgBar
  ))
  
}