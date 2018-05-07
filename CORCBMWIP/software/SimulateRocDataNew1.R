#
# simulate correlated data; this works for both 1J and IJ cases; critical contribution by XZ
# TEMP*DPC checked carefully that all values are correct

SimulateRocDataNew <- function(I, J, K1, K2, ChiNew, DesiredNumBins) {
  
  # replace missing values with values from ChiNewBar with appropriate p
  # this is the alpha_sub_t vector
  A <- array(dim= c(2, I * J))
  for (t in 1:2){
    count <- 1
    for (i in 1:I){
      for (j in 1:J){
        if (t == 1){
          A[t, count] <- 0
        }else{
          A[t, count] <- ChiNew[i, i, j, j, 1] # using diagonal elements only
        }
        count <- count + 1
      }
    }
  }
  
  SigmaA <- array(dim = c(2, I, J, I, J))
  for(t in 1:2){
    for (i1 in 1:I) {
      for (j1 in 1:J) {
        for (i2 in 1:I) {        
          for (j2 in 1:J) {
            if (t == 1){
              SigmaA[t,i1,j1,i2,j2] <- ChiNew[i1, i2, j1, j2, 2] * ChiNew[i1, i2, j1, j2, 4] * ChiNew[i1, i2, j1, j2, 5]
            }else{
              SigmaA[t,i1,j1,i2,j2] <- 1 * 1 * ChiNew[i1, i2, j1, j2, 6]
            }
          }
        }
      }
    }
  }
  
#   for (t in 1:2){
#     for (i1 in 1:I) {
#       for (j1 in 1:J) {
#         for (i2 in 1:I) {        
#           for (j2 in 1:J) {
#             cat(t, i1, j1, i2, j2,SigmaA[t,i1,j1,i2,j2],"\n")
#           }
#         }
#       }
#     }
#   }
  
  
  K <- c(K1,K2)
  tempZ <- as.list(rep(NA, 2))
  for (t in 1:2){
    tempZ[[t]] <- mvrnorm(K[t], A[t, ], nearPD(SigmaA[t, , ]))
  }
  
  z1b <- array(dim = c(I, J, K1))
  z2b <- array(dim = c(I, J, K2))
  
  count <- 1
  for (i in 1:I){
    for (j in 1:J){
      zb <- ToIntegerRatings(tempZ[[1]][, count], tempZ[[2]][, count], DesiredNumBins)
      z1b[i, j, ] <- zb$f
      z2b[i, j, ] <- zb$t
      count <- count + 1
    }
  }
  
  if (any(is.na(z1b))) stop("error z1b")
  if (any(is.na(z2b))) stop("error z2b")
  
  done <- 1
  return(list(
    z1b = z1b,
    z2b = z2b
  ))
}
