SigmaChiNew <- function(MyFile, S, I, J, ChiOrgBar, SigmaChiOrgBar)
{
  ChiNew <- array(dim = c(S, I, I, J, J, 6))
  for (i1 in 1:I) {
    for (i2 in 1:I) {        
      for (j1 in 1:J) {
        for (j2 in 1:J) {
          p <- PairingType(i1, i2, j1, j2)
          if (resampleChi){
            ChiNew[, i1, i2, j1, j2, ] <- mvrnorm(S, ChiOrgBar[p, ], Sigma = nearPD(SigmaChiOrgBar[p, , ]))
          }else{
            ChiNew[, i1, i2, j1, j2, ] <- matrix(ChiOrgBar[p, ], nrow = S, ncol = 6, byrow = TRUE)
          }
        }
      }
    }
  }
  return(ChiNew)
}
